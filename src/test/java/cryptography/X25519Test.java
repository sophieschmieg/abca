package cryptography;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.IOException;

import org.junit.jupiter.api.Test;

import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.integers.Integers;
import fields.integers.Integers.IntE;

class X25519Test {

	private byte[] decodeHex(String hex) throws IOException {
		if (hex.length() % 2 != 0) {
			throw new IOException("Not an even number of chars");
		}
		byte[] result = new byte[hex.length() / 2];
		for (int i = 0; i < result.length; i++) {
			result[i] = (byte) Integer.parseInt(hex.substring(2 * i, 2 * i + 2), 16);
		}
		return result;
	}

	private IntE readHex(String hex) throws IOException {
		if (hex.length() != 64) {
			throw new IOException("wrong size");
		}
		Integers z = Integers.z();
		IntE result = z.zero();
		for (int i = 0; i < 32; i++) {
			result = z.add(result, z.multiply(z.getInteger(Integer.decode("0x" + hex.substring(2 * i, 2 * i + 2))),
					z.power(z.getInteger(256), i)));
		}
		return result;
	}

	private boolean testXCoordinate(IntE x) {
		Integers z = Integers.z();
		IntE prime = z.subtract(z.power(z.getInteger(2), 255), z.getInteger(19));
		PrimeField fp = PrimeField.getPrimeField(prime);
		PFE a = fp.getInteger(486662);
		PFE xr = fp.reduce(x);
		PFE rhs = fp.add(fp.power(xr, 3), fp.multiply(a, xr, xr), xr);
		return fp.hasSqrt(rhs);
	}

	@Test
	void testTwist() throws IOException {
		Integers z = Integers.z();
		IntE x = z.parse("88838573511839298940907593866106493194" + "17338800022198945255395922347792736741");
		assertFalse(testXCoordinate(x));
		String hex = "eaffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff";
		IntE r = readHex(hex);
		assertTrue(testXCoordinate(r));
		String hex2 = "eaffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff7f";
		IntE r2 = readHex(hex2);
		assertFalse(testXCoordinate(r2));
		String hex3 = "dbffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff";
		IntE r3 = readHex(hex3);
		assertTrue(testXCoordinate(r3));
		String hex4 = "dbffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff7f";
		IntE r4 = readHex(hex4);
		assertFalse(testXCoordinate(r4));
		String hex5 = "dcffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff";
		IntE r5 = readHex(hex5);
		assertFalse(testXCoordinate(r5));
		String hex6 = "dcffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff7f";
		IntE r6 = readHex(hex6);
		assertTrue(testXCoordinate(r6));
		String hex7 = "daffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff";
		IntE r7 = readHex(hex7);
		assertTrue(testXCoordinate(r7));
		String hex8 = "daffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff7f";
		IntE r8 = readHex(hex8);
		assertTrue(testXCoordinate(r8));
		String hex9 = "e5210f12786811d3f4b7959d0538ae2c31dbe7106fc03c3efc4cd549c715a493";
		IntE r9 = readHex(hex9);
		assertTrue(testXCoordinate(r9));
		String hex10 = "e5210f12786811d3f4b7959d0538ae2c31dbe7106fc03c3efc4cd549c715a413";
		IntE r10 = readHex(hex10);
		assertFalse(testXCoordinate(r10));
	}

	@Test
	void test() throws IOException {
		X25519 x25519 = new X25519(decodeHex("a546e36bf0527c9d3b16154b82465edd62144c0ac1fc5a18506a2244ba449ac4"));
		ByteArray result = x25519
				.evaluate(new ByteArray(decodeHex("e6db6867583030db3594c1a424b15f7c726624ec26b3353b10a903a6d0ab1c4c")));
		assertEquals("c3 da 55 37 9d e9 c6 90 8e 94 ea 4d f2 8d 08 4f 32 ec cf 03 49 1c 71 f7 54 b4 07 55 77 a2 85 52",
				result.toString());
		x25519 = new X25519(decodeHex("4b66e9d4d1b4673c5ad22691957d6af5c11b6421e0ea01d42ca4169e7918ba0d"));
		result = x25519
				.evaluate(new ByteArray(decodeHex("e5210f12786811d3f4b7959d0538ae2c31dbe7106fc03c3efc4cd549c715a493")));
		assertEquals("95 cb de 94 76 e8 90 7d 7a ad e4 5c b4 b8 73 f8 8b 59 5a 68 79 9f a1 52 e6 f8 f7 64 7a ac 79 57",
				result.toString());
	}
}
