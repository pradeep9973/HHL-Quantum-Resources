"prints a list of eigenvalues and their inverted counterpart(101 and 010) for GGM calculation"


numbers = []
for i in range(1, 101):
    binary = bin(i)[2:]
    flipped = ''
    for bit in binary:
        if bit == '0':
            flipped += '1'
        else:
            flipped += '0'
    decimal_flipped = int(flipped, 2)
    if decimal_flipped != 0:
        ratio = i/decimal_flipped
        numbers.append((i, decimal_flipped, ratio, binary, flipped))

numbers.sort(key=lambda x: x[2])

print("Number".ljust(10) + "Flipped Decimal".ljust(15) + "Ratio".ljust(10) + "Binary".ljust(15) + "Flipped Binary".ljust(15))
for number in numbers:
    print(str(number[0]).ljust(10) + str(number[1]).ljust(15) + format(number[2], '.3f').ljust(10) + number[3].ljust(15) + number[4].ljust(15))
