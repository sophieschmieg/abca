package cryptography;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;

import util.Pair;

class KyberTest {

	@Test
	void kyberTest() {
		Kyber kyber = new Kyber(2, 3, 2, 10, 4, Sha3.SHA3_512, Sha3.SHA3_256, Sha3.SHAKE_256_PRF, Sha3.SHAKE_128,
				Sha3.SHAKE_256);
		ByteArray privateKey = kyber.createPrivateKey();
		ByteArray publicKey = kyber.createPublicKey(privateKey);
		System.out.println(privateKey);
		System.out.println(publicKey);
		Pair<VariableLengthKey, ByteArray> encapsulated = kyber.encapsulate(publicKey);
		System.out.println(encapsulated.getSecond());
		byte[] aliceKey = encapsulated.getFirst().key(32);
		System.out.println(new ByteArray(aliceKey));
		VariableLengthKey decrypted = kyber.decapsulate(encapsulated.getSecond(), privateKey);
		byte[] bobKey = decrypted.key(32);
		System.out.println(new ByteArray(bobKey));
		assertEquals(new ByteArray(aliceKey), new ByteArray(bobKey));
	}

	@Test
	void funSizeKyberTest() {
		Kyber kyber = new Kyber(8, 41, 3, 2, 1, 1, 6, 6, Sha2.SHA2_512, Sha2.SHA2_256, new Hmac(Sha2.SHA2_512), new Hkdf(new Hmac(Sha2.SHA2_512)),
				new Hkdf(new Hmac(Sha2.SHA2_512), new byte[0x01]));
		ByteArray privateKey = kyber.createPrivateKey();
		ByteArray publicKey = kyber.createPublicKey(privateKey);
		System.out.println(privateKey);
		System.out.println(publicKey);
		Pair<VariableLengthKey, ByteArray> encapsulated = kyber.encapsulate(publicKey);
		System.out.println(encapsulated.getSecond());
		byte[] aliceKey = encapsulated.getFirst().key(32);
		System.out.println(new ByteArray(aliceKey));
		VariableLengthKey decrypted = kyber.decapsulate(encapsulated.getSecond(), privateKey);
		byte[] bobKey = decrypted.key(32);
		System.out.println(new ByteArray(bobKey));
		assertEquals(new ByteArray(aliceKey), new ByteArray(bobKey));
	}

}
