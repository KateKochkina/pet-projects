import sys

digit_string = sys.argv[1]
sum = 0

for digit in digit_string:
    sum += int(digit)

print(sum)
