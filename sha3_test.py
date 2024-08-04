import binascii
import hashlib
import time

def test_hexdigit(ch):
    if '0' <= ch <= '9':
        return ord(ch) - ord('0')
    if 'A' <= ch <= 'F':
        return ord(ch) - ord('A') + 10
    if 'a' <= ch <= 'f':
        return ord(ch) - ord('a') + 10
    return -1

def test_readhex(hex_str, maxbytes):
    """
    Convert a hexadecimal string into a byte array.
    1. Read a hexadecimal string.
    2. Convert it into a sequence of bytes.
    3. Return the resulting byte array.
    """
    buf = bytearray()
    for i in range(min(maxbytes, len(hex_str) // 2)):
        h = test_hexdigit(hex_str[2 * i])
        if h < 0:
            break
        l = test_hexdigit(hex_str[2 * i + 1])
        if l < 0:
            break
        buf.append((h << 4) + l)
    return buf

def test_sha3():
    testvec = [
        ("", "6B4E03423667DBB73B6E15454F0EB1ABD4597F9A1B078E3F5B5A6BC7"),
        ("9F2FCC7C90DE090D6B87CD7E9718C1EA6CB21118FC2D5DE9F97E5DB6AC1E9C10",
         "2F1A5F7159E34EA19CDDC70EBF9B81F1A66DB40615D7EAD3CC1F1B954D82A3AF"),
        ("E35780EB9799AD4C77535D4DDB683CF33EF367715327CF4C4A58ED9CBDCDD486F669F80189D549A9364FA82A51A52654EC721BB3AAB95DCEB4A86A6AFA93826DB923517E928F33E3FBA850D45660EF83B9876ACCAFA2A9987A254B137C6E140A21691E1069413848",
         "D1C0FA85C8D183BEFF99AD9D752B263E286B477F79F0710B010317017397813344B99DAF3BB7B1BC5E8D722BAC85943A"),
        ("3A3A819C48EFDE2AD914FBF00E18AB6BC4F14513AB27D0C178A188B61431E7F5623CB66B23346775D386B50E982C493ADBBFC54B9A3CD383382336A1A0B2150A15358F336D03AE18F666C7573D55C4FD181C29E6CCFDE63EA35F0ADF5885CFC0A3D84A2B2E4DD24496DB789E663170CEF74798AA1BBCD4574EA0BBA40489D764B2F83AADC66B148B4A0CD95246C127D5871C4F11418690A5DDF01246A0C80A43C70088B6183639DCFDA4125BD113A8F49EE23ED306FAAC576C3FB0C1E256671D817FC2534A52F5B439F72E424DE376F4C565CCA82307DD9EF76DA5B7C4EB7E085172E328807C02D011FFBF33785378D79DC266F6A5BE6BB0E4A92ECEEBAEB1",
         "6E8B8BD195BDD560689AF2348BDC74AB7CD05ED8B9A57711E9BE71E9726FDA4591FEE12205EDACAF82FFBBAF16DFF9E702A708862080166C2FF6BA379BC7FFC2")
    ]

    fails = 0
    for i, (msg_str, sha_str) in enumerate(testvec):
        msg = test_readhex(msg_str, 256)
        sha = test_readhex(sha_str, 64)

        if len(sha) == 28:
            digest = hashlib.sha3_224(msg).digest()
        elif len(sha) == 32:
            digest = hashlib.sha3_256(msg).digest()
        elif len(sha) == 48:
            digest = hashlib.sha3_384(msg).digest()
        elif len(sha) == 64:
            digest = hashlib.sha3_512(msg).digest()

        if digest != sha:
            print(f"[{i}] SHA3-{len(sha) * 8}, len {len(msg)} test FAILED.")
            fails += 1
        print(f"[{i}] {digest}, {msg}")

    return fails

def test_shake():
    testhex = [
        "43E41B45A653F2A5C4492C1ADD544512DDA2529833462B71A41A45BE97290B6F",
        "AB0BAE316339894304E35877B0C28A9B1FD166C796B9CC258A064A8F57E27F2A",
        "44C9FB359FD56AC0A9A75A743CFF6862F17D7259AB075216C0699511643B6439",
        "6A1A9D7846436E4DCA5728B6F760EEF0CA92BF0BE5615E96959D767197A0BEEB"
    ]

    fails = 0

    for i, hex_str in enumerate(testhex):
        ref = test_readhex(hex_str, 32)

        if i % 2 == 0:
            shake = hashlib.shake_128()
        else:
            shake = hashlib.shake_256()

        if i >= 2:
            pattern = bytearray([0xA3] * 200)
            shake.update(pattern)

        shake_out = shake.digest(512)
        buf = shake_out[480:512]

        if buf != ref:
            print(f"[{i}] SHAKE{128 if i % 2 == 0 else 256}, len {1600 if i >= 2 else 0} test FAILED.")
            fails += 1

    return fails

def test_speed():
    start_time = time.time()
    n = 0
    while (time.time() - start_time) < 3:
        for _ in range(100000):
            hashlib.shake_128().digest(256)
        n += 100000

    elapsed = time.time() - start_time
    keccak_per_sec = n / elapsed

    print(f"Keccak-p[1600,24] / Second: {keccak_per_sec:.3f}")

def main():
    if test_sha3() == 0 and test_shake() == 0:
        print("FIPS 202 / SHA3, SHAKE128, SHAKE256 Self-Tests OK!")
    test_speed()

if __name__ == "__main__":
    main()
