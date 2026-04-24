cd /volatile/Git/GeoLab/build

echo "Starting verification of all apps..." > ../verify_apps.log

for f in ../apps/*.cxx; do
    sourceName=$(basename "$f" .cxx)
    echo "Building $sourceName..."
    if ! make "$sourceName" -j8 > /tmp/build_${sourceName}.log 2>&1; then
        echo "FAILED: $sourceName" >> ../verify_apps.log
        grep -i "error" /tmp/build_${sourceName}.log | grep -v "Erreur 2" >> ../verify_apps.log
    else
        echo "SUCCESS: $sourceName" >> ../verify_apps.log
    fi
done

cat ../verify_apps.log
