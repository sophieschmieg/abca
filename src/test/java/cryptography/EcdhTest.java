package cryptography;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.math.BigInteger;

import org.junit.jupiter.api.Test;

import cryptography.interfaces.DhScheme.Role;
import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.integers.Integers.IntE;
import varieties.curves.elliptic.EllipticCurve;
import varieties.projective.ProjectivePoint;

class EcdhTest {
	private PrimeField fp256;
	private PFE p256A;
	private PFE p256B;
	private EllipticCurve<PFE> p256;
	private ProjectivePoint<PFE> p256Generator;
	private BigInteger p256Order;

	@Test
	void test() {
		BigInteger p256Prime = new BigInteger("ffffffff00000001000000000000000000000000ffffffffffffffffffffffff", 16);
		fp256 = PrimeField.getPrimeField(p256Prime);
		p256A = fp256.getElement(new BigInteger("-3"));
		p256B = fp256.getElement(
				new BigInteger("41058363725152142129326129780047268409114441015993725554835256314039467401291"));
		PFE baseX = fp256.getElement(
				new BigInteger("48439561293906451759052585252797914202762949526041747995844080717082404635286"));
		PFE baseY = fp256.getElement(
				new BigInteger("36134250956749795798585127919587881956611106672985015071877198253568414405109"));
		p256Generator = new ProjectivePoint<>(fp256, baseX, baseY, fp256.one());
		p256Order = new BigInteger("115792089210356248762697446949407573529996955224135760342422259061068512044369");
		p256 = new EllipticCurve<>(fp256, p256A, p256B);
		p256.setNumberOfElements(p256Order);
		Ecdh<PFE> ecdh = new Ecdh<>(fp256, p256, p256Generator);
		IntE alicePrivateKey = ecdh.createPrivateKey(Role.Alice);
		System.out.println("Private key Alice: " + alicePrivateKey);
		ProjectivePoint<PFE> alicePublicKey = ecdh.createPublicKey(alicePrivateKey);
		System.out.println("Public key Alice: " + alicePublicKey);
		IntE bobPrivateKey = ecdh.createPrivateKey(Role.Bob);
		System.out.println("Private key Bob: " + bobPrivateKey);
		ProjectivePoint<PFE> bobPublicKey = ecdh.createPublicKey(bobPrivateKey);
		System.out.println("Public key Bob: " + bobPublicKey);

		ProjectivePoint<PFE> sharedSecretAlice = ecdh.agree(alicePrivateKey, bobPublicKey);
		ProjectivePoint<PFE> sharedSecretBob = ecdh.agree(bobPrivateKey, alicePublicKey);

		System.out.println("shared secret Alice: " + sharedSecretAlice);
		System.out.println("shared secret Bob:   " + sharedSecretBob);
		assertEquals(sharedSecretAlice, sharedSecretBob);
	}

}
