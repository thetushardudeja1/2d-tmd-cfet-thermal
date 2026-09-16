"""Command runner for a compute host reached through an SSH jump host (ProxyJump).

The target compute host is LAN-only (no public route), so this hops through a
public-facing jump host with a direct-tcpip tunnel, then authenticates to the
target over that tunnel.

    python ssh_run_atq.py "<command>"

Examples:
    python ssh_run_atq.py "hostname; uptime; nproc"
    python ssh_run_atq.py "source ~/.bashrc; which atkpython"

Notes:
  * Non-interactive bash by default. The target's software environment may
    live in ~/.cshrc; if a command needs the csh environment, wrap it:
    csh -c 'source ~/.cshrc; ...'
  * All connection details (hosts, user, passwords) are read from ./config.txt
    (gitignored, see config.example.txt for the format). Keep config.txt in
    this project folder, beside this script. Nothing here is hardcoded.
  * Short 60 s command timeout. For multi-hour jobs, launch detached
    (setsid nohup ... &) and poll the log file instead of waiting on this call.
"""
import os
import sys
import paramiko

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
sys.stderr.reconfigure(encoding="utf-8", errors="replace")

# --- filled in by load_config() from config.txt (gitignored, not in this repo) ---
JUMP_HOST = None   # public entry / jump host
JUMP_USER = None
JUMP_PASS = None

TGT_HOST = None    # LAN-only target compute host
TGT_PORT = 22
TGT_USER = None
TGT_PASS = None    # config.txt "atq:" line


def load_config():
    """Read host/user/pwd/atq/tgt_host from config.txt beside this script."""
    global JUMP_HOST, JUMP_USER, JUMP_PASS, TGT_HOST, TGT_USER, TGT_PASS
    cfg = os.path.join(os.path.dirname(os.path.abspath(__file__)), "config.txt")
    if not os.path.exists(cfg):
        raise SystemExit(
            "config.txt not found. Copy config.example.txt to config.txt and "
            "fill in host/user/pwd/atq/tgt_host (this file is gitignored, never commit it)."
        )
    vals = {}
    with open(cfg, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if ":" in line:
                k, v = line.split(":", 1)
                vals[k.strip().lower()] = v.strip()
    JUMP_HOST = vals.get("host", JUMP_HOST)
    JUMP_USER = vals.get("user", JUMP_USER)
    JUMP_PASS = vals.get("pwd", JUMP_PASS)
    TGT_HOST = vals["tgt_host"]
    TGT_USER = vals.get("user", JUMP_USER)
    TGT_PASS = vals.get("atq", TGT_PASS)


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    cmd = sys.argv[1]
    load_config()

    jump = paramiko.SSHClient()
    jump.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    jump.connect(JUMP_HOST, username=JUMP_USER, password=JUMP_PASS,
                 timeout=25, allow_agent=False, look_for_keys=False)

    chan = jump.get_transport().open_channel(
        "direct-tcpip", (TGT_HOST, TGT_PORT), ("127.0.0.1", 0), timeout=25)

    tgt = paramiko.SSHClient()
    tgt.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    tgt.connect(TGT_HOST, port=TGT_PORT, username=TGT_USER, password=TGT_PASS,
                sock=chan, timeout=25, allow_agent=False, look_for_keys=False)

    _, out, err = tgt.exec_command(cmd, timeout=60)
    o = out.read().decode(errors="replace")
    e = err.read().decode(errors="replace")
    print(o)
    if e:
        print("STDERR:", e, file=sys.stderr)

    tgt.close()
    jump.close()


if __name__ == "__main__":
    main()
