#!/usr/bin/env bash
ArrayExpress_location="$@";
grep 'transcription profiling by array' "$@" | awk '{print $1}'
