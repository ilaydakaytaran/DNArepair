import sys

def main():
  from cbio_py import cbio_mod as cb 
  d = cb.getAliasForGene(2705)
  print(d)


if __name__ == "__main__":
  if len(sys.argv) != 1:
    print("Usage failed: python trial.py")
    sys.exit(1)
  main()

