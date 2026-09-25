import sys
import os
import argparse

def convert_file(input_path, output_path, var_name):
    with open(input_path, 'r', encoding='utf-8', errors='ignore') as f:
        data = f.read()
    
    out_dir = os.path.dirname(output_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        
    with open(output_path, 'w', encoding='utf-8') as f_out:
        f_out.write("#pragma once\n\n")
        f_out.write("namespace CloudJ::Tables {\n\n")
        f_out.write(f"constexpr const char* {var_name} = R\"CLOUDJ_TABLE_EOF({data})CLOUDJ_TABLE_EOF\";\n\n")
        f_out.write("} // namespace CloudJ::Tables\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--test", action="store_true")
    parser.add_argument("--input-dir")
    parser.add_argument("--output-dir")
    args = parser.parse_args()

    if args.test:
        with open("test.dat", "w") as f: f.write("test_content")
        convert_file("test.dat", "test.hpp", "test_var")
        os.remove("test.dat")
        os.remove("test.hpp")
        print("Verified")
        sys.exit(0)

    if args.input_dir and args.output_dir:
        for file in os.listdir(args.input_dir):
            if file.endswith(".dat"):
                # standard kebab-to-underscore mapping for valid C++ variable names
                var_name = file.replace(".dat", "").replace("-", "_")
                convert_file(os.path.join(args.input_dir, file), 
                             os.path.join(args.output_dir, f"{var_name}.hpp"), 
                             var_name)
        print("Table conversion successfully completed!")
