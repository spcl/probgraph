import os
import argparse
import math

Header = "type,Problem,approximation-scheme,baseline,graph-name,thread-count,vertices,edges,treshold,b,m,preprocessing-time,tc-time,total-runtime,size-of-BF,CSR-size,total-size,approximated-count\n"

def make_csv(directory, outfilename = "results.csv"):
    with open(outfilename, 'w') as outfile:
        outfile.writelines(Header)
        
        for subdir, dirs, files in os.walk(directory):
            for file in files:
                filename = os.path.join(subdir, file);
                if "gather" in filename:
                    with open(filename, 'r') as infile:
                        
                        line = " ";                        
                        while ( line ):
                            line = infile.readline();
                            if  "RRR" in line:
                                tmp_RRR_line = line.split(" ")[0:-1]

                                #hotfix for jp output problem 
                                if tmp_RRR_line[1].split("-")[0] == "JP":
                                    RRR_line = tmp_RRR_line[0:11]
                                    RRR_line.extend(tmp_RRR_line[12:])
                                else:
                                    RRR_line = tmp_RRR_line


                            if  "SSS" in line:
                                SSS_line = line.split(" ")[-4:]
                        merged_line = ",".join(RRR_line + SSS_line)
                        outfile.writelines(merged_line);


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Collect results for approximate algorithms into a csv")
    
    parser.add_argument("--input-dir", default="/",
        help="directory that contains the already run experiments")
    parser.add_argument("--output-name", default="/",
        help="output file")

    args = parser.parse_args()

    input_dir = args.input_dir;
    csv_filename = args.output_name;
     
    make_csv(input_dir, csv_filename);


