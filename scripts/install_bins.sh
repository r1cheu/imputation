#!/usr/bin/env bash
# Install native bioinformatics binaries to <project>/bin without conda.
#
# Portability:
#   - bwa-mem2 / fastp: official prebuilt static binaries.
#   - htslib / bcftools / samtools: built from source. htslib has runtime
#     CPU dispatch (__builtin_cpu_supports + __attribute__((target))), so
#     SIMD differences across nodes don't matter. The only portability axis
#     is glibc -- run this script on the OLDEST node you target so the
#     resulting binaries link against an old enough glibc.
#
# Build cache: /tmp/$USER/imputation-bins (local disk, not NFS).
#
# Usage:
#   bash scripts/install_bins.sh                # install all
#   bash scripts/install_bins.sh bwa-mem2 fastp # install selected
#   FORCE=1 bash scripts/install_bins.sh        # rebuild even if installed
#
# NOTE: GLIMPSE1 is not installed here -- it is pre-built on the cluster
# and expected to be on PATH (GLIMPSE_chunk / GLIMPSE_phase / GLIMPSE_ligate).

set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
BIN="$ROOT/bin"
CACHE="/tmp/${USER:-build}/imputation-bins"
JOBS=${JOBS:-$(nproc 2>/dev/null || echo 4)}
FORCE=${FORCE:-0}

BWA_MEM2_VER=2.3
FASTP_VER=0.23.4
HTSLIB_VER=1.23.1
BCFTOOLS_VER=1.23.1
SAMTOOLS_VER=1.23.1

mkdir -p "$BIN" "$CACHE"

c_reset=$'\033[0m'; c_cyan=$'\033[36m'; c_green=$'\033[32m'; c_red=$'\033[31m'; c_yellow=$'\033[33m'
log()  { printf '%s[install-bins]%s %s\n' "$c_cyan"   "$c_reset" "$*"; }
ok()   { printf '%s[ok]%s %s\n'           "$c_green"  "$c_reset" "$*"; }
warn() { printf '%s[warn]%s %s\n'         "$c_yellow" "$c_reset" "$*" >&2; }
die()  { printf '%s[error]%s %s\n'        "$c_red"    "$c_reset" "$*" >&2; exit 1; }

need_install() {
    local target=$1
    if [[ "$FORCE" == 1 ]]; then return 0; fi
    [[ ! -x "$BIN/$target" ]]
}

fetch() {
    local url=$1 out=$2
    if [[ -s "$out" ]]; then return; fi
    log "downloading $(basename "$out")"
    curl -fL --no-progress-meter --retry 3 --retry-delay 2 -o "$out.part" "$url"
    mv "$out.part" "$out"
}

# ---------- 1. bwa-mem2 (prebuilt with runtime SIMD dispatch) ----------
install_bwa_mem2() {
    if ! need_install bwa-mem2; then ok "bwa-mem2 $BWA_MEM2_VER already installed"; return; fi
    local ver=$BWA_MEM2_VER
    local tarball="bwa-mem2-${ver}_x64-linux.tar.bz2"
    local url="https://github.com/bwa-mem2/bwa-mem2/releases/download/v${ver}/${tarball}"
    local src="$CACHE/$tarball"
    fetch "$url" "$src"

    log "extracting bwa-mem2 $ver"
    local d="$CACHE/bwa-mem2-${ver}"
    rm -rf "$d"
    mkdir -p "$d"
    tar -xjf "$src" -C "$d" --strip-components=1
    install -m755 "$d/bwa-mem2" "$BIN/bwa-mem2"
    for variant in avx512bw avx2 avx sse42 sse41; do
        if [[ -f "$d/bwa-mem2.$variant" ]]; then
            install -m755 "$d/bwa-mem2.$variant" "$BIN/bwa-mem2.$variant"
        fi
    done
    ok "bwa-mem2 $ver"
}

# ---------- 2. fastp (fully static prebuilt) ----------
install_fastp() {
    if ! need_install fastp; then ok "fastp $FASTP_VER already installed"; return; fi
    local ver=$FASTP_VER
    local url="http://opengene.org/fastp/fastp.${ver}"
    local src="$CACHE/fastp-$ver"
    fetch "$url" "$src"
    install -m755 "$src" "$BIN/fastp"
    ok "fastp $ver"
}

# ---------- htslib (built once, shared by samtools+bcftools) ----------
# Returns the htslib build dir path; configure --with-htslib=<dir> wants a
# built source tree containing libhts.a + htslib/ headers.
build_htslib() {
    local ver=$HTSLIB_VER
    local d="$CACHE/htslib-${ver}"
    if [[ -f "$d/libhts.a" && -x "$d/tabix" ]]; then
        echo "$d"; return
    fi
    fetch "https://github.com/samtools/htslib/releases/download/${ver}/htslib-${ver}.tar.bz2" "$CACHE/htslib-${ver}.tar.bz2"

    log "building htslib $ver"
    rm -rf "$d"
    mkdir -p "$d"
    tar -xjf "$CACHE/htslib-${ver}.tar.bz2" -C "$d" --strip-components=1
    (
        cd "$d"
        ./configure --disable-libcurl --disable-s3 --disable-gcs >/dev/null
        make -j"$JOBS" >/dev/null
    )
    echo "$d"
}

# ---------- 3. bcftools ----------
install_bcftools() {
    if ! need_install bcftools; then ok "bcftools $BCFTOOLS_VER already installed"; return; fi
    local ver=$BCFTOOLS_VER
    local hts; hts=$(build_htslib)
    fetch "https://github.com/samtools/bcftools/releases/download/${ver}/bcftools-${ver}.tar.bz2" "$CACHE/bcftools-${ver}.tar.bz2"

    log "building bcftools $ver"
    local d="$CACHE/bcftools-${ver}"
    rm -rf "$d"
    mkdir -p "$d"
    tar -xjf "$CACHE/bcftools-${ver}.tar.bz2" -C "$d" --strip-components=1
    (
        cd "$d"
        ./configure --with-htslib="$hts" >/dev/null
        make -j"$JOBS" >/dev/null
    )
    install -m755 "$d/bcftools" "$BIN/bcftools"
    install -m755 "$hts/tabix"  "$BIN/tabix"
    install -m755 "$hts/bgzip"  "$BIN/bgzip"
    ok "bcftools $ver (+ tabix, bgzip)"
}

# ---------- 4. samtools ----------
install_samtools() {
    if ! need_install samtools; then ok "samtools $SAMTOOLS_VER already installed"; return; fi
    local ver=$SAMTOOLS_VER
    local hts; hts=$(build_htslib)
    fetch "https://github.com/samtools/samtools/releases/download/${ver}/samtools-${ver}.tar.bz2" "$CACHE/samtools-${ver}.tar.bz2"

    log "building samtools $ver"
    local d="$CACHE/samtools-${ver}"
    rm -rf "$d"
    mkdir -p "$d"
    tar -xjf "$CACHE/samtools-${ver}.tar.bz2" -C "$d" --strip-components=1
    (
        cd "$d"
        ./configure --with-htslib="$hts" --without-curses >/dev/null
        make -j"$JOBS" >/dev/null
    )
    install -m755 "$d/samtools" "$BIN/samtools"
    ok "samtools $ver"
}

# ---------- driver ----------
declare -A INSTALLERS=(
    [bwa-mem2]=install_bwa_mem2
    [fastp]=install_fastp
    [bcftools]=install_bcftools
    [samtools]=install_samtools
)

ORDER=(bwa-mem2 fastp bcftools samtools)

if (( $# == 0 )); then
    targets=("${ORDER[@]}")
else
    targets=("$@")
fi

for t in "${targets[@]}"; do
    fn=${INSTALLERS[$t]:-}
    [[ -z "$fn" ]] && die "unknown target: $t (available: ${!INSTALLERS[*]})"
    log "==> $t"
    "$fn"
done

log "done. binaries in $BIN"
log "add to PATH:  export PATH=\"$BIN:\$PATH\""
