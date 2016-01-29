# Convert .tsv file to .csv
# Run in /bin directory to work

tr '\t' ',' < ../raw-data/CosmicFusionExport_v69_310514.tsv > \
    ../data/CosmicFusionExport_v69_310514.csv
