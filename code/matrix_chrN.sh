#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_matrix_file> <output_matrix_file>"
    echo "Example: $0 ../ctrl1.matrix ctrl1_converted.matrix"
    exit 1
fi

INPUT_FILE="$1"
OUTPUT_FILE="$2"

if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file '$INPUT_FILE' not found."
    exit 1
fi

# Create the EXTENDED mapping file
MAPPING_FILE=$(mktemp)

cat << 'EOF' > "$MAPPING_FILE"
chr1	NC_000067.7
chr2	NC_000068.8
chr3	NC_000069.7
chr4	NC_000070.7
chr5	NC_000071.7
chr6	NC_000072.7
chr7	NC_000073.7
chr8	NC_000074.7
chr9	NC_000075.7
chr10	NC_000076.7
chr11	NC_000077.7
chr12	NC_000078.7
chr13	NC_000079.7
chr14	NC_000080.7
chr15	NC_000081.7
chr16	NC_000082.7
chr17	NC_000083.7
chr18	NC_000084.7
chr19	NC_000085.7
chrX	NC_000086.8
chrY	NC_000087.8
chrM	NC_005089.1
chrUn_NT_162750	NT_162750.1
chrUn_NT_165789	NT_165789.3
chrUn_NT_166280	NT_166280.1
chrUn_NT_166281	NT_166281.1
chrUn_NT_166282	NT_166282.1
chrUn_NT_166307	NT_166307.1
chrUn_NT_166438	NT_166438.1
chrUn_NT_166443	NT_166443.1
chrUn_NT_166450	NT_166450.1
chrUn_NT_166456	NT_166456.1
chrUn_NT_166462	NT_166462.1
chrUn_NT_187055	NT_187055.1
chrUn_NT_187056	NT_187056.1
chrUn_NT_187057	NT_187057.1
chrUn_NT_187058	NT_187058.1
chrUn_NT_187059	NT_187059.1
chrUn_NT_187060	NT_187060.1
chrUn_NT_187064	NT_187064.1
chrUn_NW_023337852	NW_023337852.1
chrUn_NW_023337853	NW_023337853.1
EOF


echo "Processing Matrix: $INPUT_FILE -> $OUTPUT_FILE"

awk '
    BEGIN { OFS="\t" } 
    
    # Load mapping (NR==FNR)
    NR==FNR {
        # We strip the version from the accession in column 2 of mapping file
        split($2, a, "."); 
        accession_base = a[1];
        
        # Map["NC_000067"] = "chr1"
        map[accession_base] = $1; 
        next;
    }

    # Process Matrix File
    {
        # Accession is in Column 3
        split($3, b, ".");
        base_key = b[1];
        
        if (base_key in map) {
            $3 = map[base_key];
        } else {
            # Optional: If you want to see what failed to map, uncomment below
            # print "Warning: No map for " $3 > "/dev/stderr"
        }
        print $0;
    }
' "$MAPPING_FILE" "$INPUT_FILE" > "$OUTPUT_FILE"

rm "$MAPPING_FILE"
echo "Done."
