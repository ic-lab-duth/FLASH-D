import subprocess as sbp

files = [
  "./src/*.h",
  "./src/main.cc"
]


files = " ".join(files)
sbp.run(
  f"tar -czvf flash_d.tar.gz {files}", shell=True
)