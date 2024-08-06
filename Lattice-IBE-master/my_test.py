import subprocess

arg1 = "Hello"
arg2 = 42

# Call the C++ program with arguments
subprocess.run(["./IBE", arg1, str(arg2)])
