import argparse
import sys

def version():
    print "0.0.1"
    sys.exit(0)

if __name__ == '__main__':
    description = "A simple test for the python gui."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--version",
                        help="Prints out the version number and exits",
                        action="store_true")
    parser.add_argument("--input",
                        help="Path to input file.",
                        default="test-data/test_in.txt")
    parser.add_argument("--output",
                        help="Full path, "
                        "including file name where output should be placed.",
                        default="test-data/test_out.txt")  
    parser.add_argument("--tool_directory",
                        help="The directory the tool currently resides in "
                        "(new in 15.03).",
                        default=None)
    parser.add_argument("--new_file_path",
                        help="onfig/galaxy.ini new_file_path value",
                        default=None)
    parser.add_argument("--tool_data_path",
                        help="config/galaxy.ini tool_data_path value",
                        default=None)
    parser.add_argument("--root_dir",
                        help="Top-level Galaxy source directory made absolute via os.path.abspath()",
                        default=None)
    parser.add_argument("--datatypes_config",
                        help="config/galaxy.ini datatypes_config value",
                        default=None)
    parser.add_argument("--user_id",
                        help="Email's numeric ID (id column of galaxy_user table in the database)",
                        default=None)
    parser.add_argument("--user_email",
                        help="User's email address",
                        default=None)
    parser.add_argument("--name",
                        help="The name value: ",
                        default=None)
    parser.add_argument("--id",
                        help="The Id value: ",
                        default=None)
    parser.add_argument("--astring", help="astring value", default=None)
    parser.add_argument("--boxst", help="boxst value", default=None)
    parser.add_argument("--smallint", help="smallint value", default=None)
    parser.add_argument("--anyint", help="anyint value", default=None)
    parser.add_argument("--smallfloat", help="smallfloat value", default=None)
    parser.add_argument("--anyfloat", help="anyfloat value", default=None)
    parser.add_argument("--truefalse", help="truefalse value", default=None)
    parser.add_argument("--christian", help="christian value", default=None)
    parser.add_argument("--user_file_type",
                        help="The declared file type value: ",
                        default=None)
    parser.add_argument("--galaxy_file_type",
                        help="Galaxy file type value: ",
                        default=None)

    args = parser.parse_args()
        

    if args.version:
        version()

    print >> sys.stderr, 'Hi, This is where you would find output to sys.error'

    print "Standard output gets displayed in the window"
    
    with open(args.input) as in_file:
        with open(args.output, "w") as out_file:
            for line in in_file:
                out_file.write(line) 

    print "tool_directory",args.tool_directory	
    print "new_file_path", args.new_file_path
    print "tool_data_path", args.tool_data_path
    print "root_dir",args.root_dir
    print "datatypes_config", args.datatypes_config
    print "user_id", args.user_id
    print "user_email", args.user_email
    print "name", args.name
    print "id", args.id
    
    print "astring", args.astring
    print "boxst", args.boxst
    print "smallint", args.smallint
    print "anyint", args.anyint
    print "smallfloat", args.smallfloat
    print "anyfloat", args.anyfloat
    print "truefalse", args.truefalse
    print "christian", args.christian

    print "user_file_type", args.user_file_type
    print "galaxy_file_type", args.galaxy_file_type

