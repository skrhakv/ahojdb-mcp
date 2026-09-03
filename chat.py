#!/usr/bin/env python3

from openai import OpenAI

# Token counter for measuring the MCP tool-result size exactly.
# gpt-5 family uses the o200k_base encoding.
try:
    import tiktoken
    _enc = tiktoken.get_encoding("o200k_base")
    def count_tokens(text):
        return len(_enc.encode(text))
except Exception:
    # Fallback if tiktoken isn't installed: rough ~4 chars/token estimate.
    def count_tokens(text):
        return len(text) // 4

SYSTEM_PROMPT = """
You are an assistant for the available MCP tools.

Use tools to answer questions.
Base factual answers on tool results.
Do not invent data or interpretations.

For questions the tools cannot answer,
politely explain which tools are available.
""".strip()

## Load config values

config = {}

for line in open("chat.conf"):
    if "=" not in line or line.strip().startswith("="):
        continue

    key, value = line.split("=", 1)
    config[key.strip()] = value.split("#", maxsplit=1)[0].strip().strip('"')

client = OpenAI(api_key=config["OPENAI_API_KEY"])
MCP_SERVER_URL = config["MCP_SERVER_URL"]

previous_response_id = None

total_tokens = 0

while True:
    user_input = input("\nYou: ")

    if user_input.lower() in ["quit", "exit"]:
        break

    response = client.responses.create(
        model="gpt-5-nano",
        instructions=SYSTEM_PROMPT,
        input=user_input,
        previous_response_id=previous_response_id,
        tools=[
            {
                "type": "mcp",
                "server_label": "mcp-workshop",
                "server_url": MCP_SERVER_URL,
                "require_approval": "never",
            }
        ],
    )

    print("Assistant:", response.output_text)

    usage = response.usage
    total_tokens += usage.total_tokens

    # Measure the exact tokens returned by the MCP server this turn.
    # Each tool invocation is an `mcp_call` item whose `output` holds the
    # raw tool result the server returned.
    mcp_response_tokens = 0
    mcp_call_count = 0
    for item in response.output:
        if getattr(item, "type", None) == "mcp_call":
            mcp_call_count += 1
            mcp_response_tokens += count_tokens(getattr(item, "output", "") or "")

    print(f"\n\n[This turn: {usage.total_tokens} tokens]")
    print(f"[Session total: {total_tokens} tokens]")
    print(
        f"[MCP response this turn: {mcp_response_tokens} tokens "
        f"across {mcp_call_count} tool call(s)]"
    )

    # Store session state.
    previous_response_id = response.id
