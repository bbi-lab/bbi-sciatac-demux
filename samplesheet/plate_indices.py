#!/usr/bin/env python3

print('== row fills first: n7 and p7 confirmed with Riza\'s plate layout')
print()
i = 0
for iplate in range(4):
  for irow in range(8):
    print('%s ' % (chr(96+irow+1)), end='')
    for icol in range(12):
      i = i + 1
      print(' %03d' % ( i ), end='')
    print()
  print() 
print()

print('== column fills first: n5 and p5 confirmed with Riza\'s plate layout')
print()
for iplate in range(4):
  for irow in range(8):
    print('%s ' % (chr(96+irow+1)), end='')
    for icol in range(12):
      print(' %03d' % (iplate*96+irow+icol*8+1), end='')
    print()
  print()
