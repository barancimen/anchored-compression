import lzma,sys
d=sys.stdin.buffer.read();f=[{'id':lzma.FILTER_LZMA2,'preset':9|lzma.PRESET_EXTREME,'lc':0,'lp':0,'pb':0}]
sys.stdout.buffer.write(lzma.compress(d,lzma.FORMAT_RAW,filters=f)if len(d)==125000 else lzma.decompress(d,lzma.FORMAT_RAW,filters=f))
