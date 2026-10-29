import builtins,sys,runpy,pathlib
original=builtins.open
root=pathlib.Path.cwd().resolve()
def traced(file, mode="r", *args, **kwargs):
 if isinstance(file,(str,pathlib.Path)) and "r" in mode and pathlib.Path(file).suffix in {".csv",".xlsx"}:
  with original(root/"verification/input_reads.txt","a") as f:f.write(str(pathlib.Path(file).resolve())+"\n")
 return original(file,mode,*args,**kwargs)
builtins.open=traced
runpy.run_path(sys.argv[1],run_name="__main__")
