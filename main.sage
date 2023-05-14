load("differential/checks.sage")
load("linear/checks.sage")

# pip install tabulate (prints data in a nice table)
from tabulate import tabulate
from tqdm import tqdm

from differential.check_pride import check_pride

from ciphers.aes import AES
from ciphers.ascon import Ascon
from ciphers.boomslang import Boomslang
from ciphers.craft import Craft
from ciphers.gift import Gift
from ciphers.iscream import iScream
from ciphers.keccak import Keccak
from ciphers.kuznechik import Kuznechik
from ciphers.led import LED
from ciphers.mantis import Mantis
from ciphers.midori import Midori
from ciphers.pride import Pride
from ciphers.prince import Prince
from ciphers.present import Present
from ciphers.rectangle import Rectangle
from ciphers.skinny import Skinny
from ciphers.streebog import Streebog
from ciphers.toy_cipher import Example1, Example2
CIPHER_LIST = [AES(), Ascon(), Boomslang(), Craft(), Gift(64), Gift(128),
               iScream(), Keccak(100), Kuznechik(), LED(), Mantis(), Midori(),
               Pride(), Prince(), Present(), Rectangle(), Skinny(64),
               Skinny(128), Streebog()]

RUN_EXPERIMENTS = {
    "differential": True,
    "linear": True
}


if RUN_EXPERIMENTS["differential"]:
    print("Checking for probability one differentials over two rounds/superboxes")
    header = ["Cipher", "Corollary 7", "Corollary 8", "Corollary 9"]
    data = []
    additional_examples = [Example1(), Example2()]
    progress = tqdm(CIPHER_LIST + additional_examples, bar_format='Progress: {bar:20}| ({n_fmt}/{total_fmt}) {postfix}')
    for cipher in progress:
        progress.set_postfix_str(f"currently checking {cipher.name}")
        c9, c10, c11 = check_corollaries(cipher)
        data.append([cipher.name, c9, c10, c11])
        if (c9, c10, c11) == (False, False, False):
            a = filter_output_differences(cipher.inverse())
            print()
            print(f"Filtered input differences for {cipher.name}")
            print(a)
            b = filter_output_differences(cipher)
            print(f"Filtered output differences for {cipher.name}")
            print(b)
    print(tabulate(data, headers=header, tablefmt="grid"))

    check_pride()

    print("Verifying examples")
    progress = tqdm(additional_examples, bar_format='Progress: {bar:20}| ({n_fmt}/{total_fmt}) {postfix}')
    for ex in progress:
        progress.set_postfix_str(f"currently checking {ex.name}")
        print()
        print(f"Differential for {ex.name} {'not ' if not ex.test_differential() else ''} verified")

if RUN_EXPERIMENTS["linear"]:
    print("Checking for perfect linear approximations over two rounds/superboxes")
    header = ["Cipher", "r=2", "r=3", "r=4"]
    data = []
    progress = tqdm(CIPHER_LIST, bar_format='Progress: {bar:20}| ({n_fmt}/{total_fmt}) {postfix}')
    for cipher in progress:
        progress.set_postfix_str(f"currently checking {cipher.name}")
        table_row = [cipher.name]
        for r in [2, 3, 4]:
            betas = linear_check(cipher, r)
            if betas == "NotChecked":
                table_row.append("/")
                tqdm.write(f"{cipher}: r={r}: Not checked")
                continue
            try:
                betas_len = len(betas)
            except OverflowError:
                betas_len = 2**betas.dimension()
            table_row.append(betas_len)
            if betas_len == 1: # beta = 0 is always a solution
                tqdm.write(f"{cipher}: r={r}: No perfect linear approximations found")
            else:
                tqdm.write(f"{cipher}: r={r}: Found {betas_len - 1} candidates for beta:")
                for b in betas[:16]: # print at most 16 values for beta
                    tqdm.write(f"-> b = 0x{ZZ(list(b)[::-1], 2).hex()}")
            data.append(table_row)
    print(tabulate(data, headers=header, tablefmt="grid"))
    print(tabulate(data, headers=header, tablefmt="latex"))
