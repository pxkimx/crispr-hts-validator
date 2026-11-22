#!/bin/bash
# Builds the download that goes on the website: the app plus "How to open.txt", in one folder.
#   macos/make_zip.sh [destination-folder]      (default: the website repository's public/downloads)
#
# Plain `ditto -c -k --keepParent`. It keeps the executable bits the .app needs, and deliberately does NOT
# pass --sequesterRsrc, which writes an AppleDouble copy of every file into a __MACOSX/ folder: 253 entries of
# junk that Windows and Linux users see when they unzip, and that the other tools' downloads do not have.
set -eu
HERE="$(cd "$(dirname "$0")/.." && pwd)"
DEST="${1:-$HOME/Documents/GitHub/paul-h-kim-site/public/downloads}"
VERSION="$(sed -n 's/^__version__ = "\(.*\)"/\1/p' "$HERE/crisprhts/__init__.py")"
APP="$HERE/CRISPR HTS Validator.app"
[ -d "$APP" ] || { echo "Build the app first: macos/build_app.sh"; exit 1; }
BUNDLED="$(/usr/bin/defaults read "$APP/Contents/Info.plist" CFBundleShortVersionString)"
[ "$BUNDLED" = "$VERSION" ] || { echo "The app says $BUNDLED but the code says $VERSION — rebuild it first."; exit 1; }

ZIP="$DEST/CRISPRHTSValidator-$VERSION-mac.zip"
STAGE="$(mktemp -d)/CRISPR HTS Validator"
mkdir -p "$STAGE"
cp -R "$APP" "$STAGE/"
cp "$HERE/macos/HOW-TO-OPEN.txt" "$STAGE/How to open.txt"
rm -f "$ZIP"
ditto -c -k --keepParent "$STAGE" "$ZIP"

if unzip -l "$ZIP" | grep -q __MACOSX; then echo "!! __MACOSX junk got in"; exit 1; fi
echo "Built: $ZIP ($(du -h "$ZIP" | cut -f1))"
echo "Remember to point the site at this file: version, size and href in lib/site-content.ts and app/page.tsx."
