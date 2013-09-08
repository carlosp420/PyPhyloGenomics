import os

def main():
    names = os.listdir(os.curdir)
    for name in names:
        if name[:5] == "test_" and name[-3:] == ".py":
            print name
            os.system("python " + name)

if __name__ == "__main__":
    main()
