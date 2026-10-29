#!/usr/bin/env bash
# download_prot_pages.sh
# Downloads every page of an InterPro API list as prot_0001.json, prot_0002.json, ...
# Each page's "next" field holds the URL of the following page; the last page has "next": null.
#
# Usage : bash download_prot_pages.sh [output_dir]      (default output_dir: prot_json)
# Needs : curl. jq is recommended (a grep/sed fallback is used when jq is missing).
# Re-running the script resumes after the last page already saved.

set -uo pipefail

START_URL='https://www.ebi.ac.uk/interpro/api/proteome/uniprot/entry/InterPro/?page_size=300&extra_fields=counters%3Aentry-protein'
OUT_DIR="${1:-prot_json}"
PREFIX="prot_"
PAUSE=1          # seconds between requests, as in InterPro's own example scripts
MAX_TRIES=10     # attempts per page before giving up

mkdir -p "$OUT_DIR" || exit 1

have_jq() { command -v jq >/dev/null 2>&1; }

# Print the "next" URL stored in a page file (prints nothing when it is null).
next_url() {
  if have_jq; then
    jq -r '.next // empty' "$1"
  else
    grep -o '"next": *"[^"]*"' "$1" | head -n 1 | sed -E 's/^"next": *"(.*)"$/\1/'
  fi
}

# Succeed only if a downloaded file is complete JSON with a "results" array.
is_complete() {
  if have_jq; then
    jq -e '.results | type == "array"' "$1" >/dev/null 2>&1
  else
    grep -q '"results"' "$1" &&
      [ "$(tail -c 64 "$1" | tr -d '[:space:]' | tail -c 1)" = "}" ]
  fi
}

# ---- resume: continue from the "next" link of the last page already saved ----
page=0
url="$START_URL"
last=""
for f in "$OUT_DIR/$PREFIX"[0-9]*.json; do
  b=$(basename "$f" .json)
  [[ ${b#"$PREFIX"} =~ ^[0-9]+$ ]] && last="$f"
done
if [ -n "$last" ]; then
  b=$(basename "$last" .json)
  page=$((10#${b#"$PREFIX"}))
  url=$(next_url "$last")
  if [ -z "$url" ]; then
    echo "Nothing to do: $last is already the last page."
    exit 0
  fi
  echo "Resuming after $last"
fi

# ---- main loop: one page per iteration until "next" is null ----
while [ -n "$url" ]; do
  page=$((page + 1))
  out=$(printf '%s/%s%04d.json' "$OUT_DIR" "$PREFIX" "$page")
  tmp="$out.part"

  try=1
  while :; do
    code=$(curl -sSL --compressed -H 'Accept: application/json' \
                --connect-timeout 30 --max-time 300 \
                -o "$tmp" -w '%{http_code}' "$url")
    if [ "$code" = "200" ] && is_complete "$tmp"; then
      break
    elif [ "$code" = "204" ]; then
      echo "HTTP 204 (no content) for page $page - stopping."
      rm -f "$tmp"
      exit 0
    elif [ "$code" = "408" ]; then
      wait_s=61    # slow query: InterPro keeps computing it in the background
    else
      wait_s=$((try * 15))
    fi
    if [ "$try" -ge "$MAX_TRIES" ]; then
      echo "Giving up on page $page (last HTTP code: $code). Re-run the script to resume." >&2
      rm -f "$tmp"
      exit 1
    fi
    echo "  page $page: HTTP $code, retry $try/$MAX_TRIES in ${wait_s}s" >&2
    sleep "$wait_s"
    try=$((try + 1))
  done

  mv "$tmp" "$out"
  if have_jq; then
    echo "Saved $out ($(jq '.results | length' "$out") records; API count $(jq '.count' "$out"))"
  else
    echo "Saved $out"
  fi
  url=$(next_url "$out")
  [ -n "$url" ] && sleep "$PAUSE"
done

echo "Finished: $page pages in $OUT_DIR/  ->  next: Rscript combine_prot_json.R $OUT_DIR"
