# define args

PARMS_FILE=$1
SCRIPT=$2
OUTPUT_FILE=$3

# Clear or create the output file
> "$OUTPUT_FILE"

# Process the input file line by line until there are no more lines
while IFS= read -r line || [ -n "$line" ]; do
    # Split the line into separate arguments and pass to Python script
    # This will preserve each space-separated value as a distinct argument
    python3 "$PYTHON_SCRIPT" $line >> "$OUTPUT_FILE"
    
    # Check if the Python script executed successfully. fi ends the if loop
    if [ $? -ne 0 ]; then
        echo "Error: Python script failed on input line: $line"
    fi
#this is a redirection operator that tells bash which file to read line by line
done < "$PARMS_FILE" 

echo "Processing complete. Results saved to $OUTPUT_FILE"