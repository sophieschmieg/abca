package main;

import java.math.BigInteger;
import java.security.KeyFactory;
import java.security.KeyPair;
import java.security.KeyPairGenerator;
import java.security.PrivateKey;
import java.security.PublicKey;
import java.security.SecureRandom;
import java.security.interfaces.ECPrivateKey;
import java.security.interfaces.ECPublicKey;
import java.security.spec.ECField;
import java.security.spec.ECFieldFp;
import java.security.spec.ECParameterSpec;
import java.security.spec.ECPoint;
import java.security.spec.ECPrivateKeySpec;
import java.security.spec.ECPublicKeySpec;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import javax.crypto.KeyAgreement;

import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import varieties.curves.EllipticCurve;
import varieties.projective.ProjectivePoint;

public class TestAsymmetricCrypto {
	private SecureRandom rand;
	private PrimeField fp256;
	private PFE p256A;
	private PFE p256B;
	private EllipticCurve<PFE> p256;
	private ProjectivePoint<PFE> p256Generator;
	private BigInteger p256Order;
	private EllipticCurve<PFE> supersingularCurve;
	private Map<Integer, ProjectivePoint<PFE>> supersingularTorsionPoint;

	public TestAsymmetricCrypto() {
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
		rand = new SecureRandom();
		supersingularTorsionPoint = new TreeMap<>();
		supersingularCurve = new EllipticCurve<>(fp256, p256A, fp256.zero());
	}

	public void doKeyAgreement() {
		BigInteger privateKeyAlice = makePrivateKey();
		ProjectivePoint<PFE> publicKeyAlice = p256.multiply(privateKeyAlice, p256Generator);
		KeyPair keyBob = makeECKeysJava();
		BigInteger privateKeyBob = ((ECPrivateKey) keyBob.getPrivate()).getS();
		ECPoint publicKeyBob = ((ECPublicKey) keyBob.getPublic()).getW();
		byte[] sharedSecretAlice = computeSharedSecret(privateKeyAlice, fromECPoint(publicKeyBob));
		byte[] sharedSecretBob = computeSharedSecretJava(privateKeyBob, toECPoint(publicKeyAlice));
		System.out.println(Arrays.equals(sharedSecretAlice, sharedSecretBob));
		for (byte b : sharedSecretAlice) {
			System.out.print(b + ", ");
		}
		System.out.println();
		for (byte b : sharedSecretBob) {
			System.out.print(b + ", ");
		}
	}

	public void doRelatedSupersingularCurveAttack(int numbits) {
		int order = 1 << numbits;
		System.out.println("Stealing " + numbits + " bits from Bob's private key:");
		KeyPair keyBob = makeECKeysJava();
		BigInteger privateKeyBob = ((ECPrivateKey) keyBob.getPrivate()).getS();
		ProjectivePoint<PFE> torsionPoint = null;

		if (supersingularTorsionPoint.containsKey(order)) {
			torsionPoint = supersingularTorsionPoint.get(order);
		} else {
			List<ProjectivePoint<PFE>> torsionPoints = supersingularCurve.getTorsionPoints(order);
			for (ProjectivePoint<PFE> p : torsionPoints) {
				if (!supersingularCurve.multiply(order / 2, p).equals(supersingularCurve.neutral())) {
					torsionPoint = p;
					supersingularTorsionPoint.put(order, torsionPoint);
					break;
				}
			}
		}
		Map<ProjectivePoint<PFE>, Integer> results = new TreeMap<>();
		for (int i = 0; i < order; i++) {
			ProjectivePoint<PFE> mult = supersingularCurve.multiply(i, torsionPoint);
			results.put(mult, i);
		}
		ProjectivePoint<PFE> sharedPoint = p256.multiply(privateKeyBob, torsionPoint);
		if (!results.containsKey(sharedPoint)) {
			System.err.println("Something went wrong!");
		} else {
			System.out.println(results.get(sharedPoint));
		}
		System.out.println(privateKeyBob.mod(BigInteger.valueOf(order)));
		System.out.println(privateKeyBob.toString(16));

	}
	
	public BigInteger makePrivateKey() {
		BigInteger privateKey;
		do {
			byte privateBytes[] = new byte[p256Order.toByteArray().length];
			rand.nextBytes(privateBytes);
			privateKey = new BigInteger(privateBytes);
		} while (privateKey.compareTo(p256Order) >= 0);
		return privateKey;
	}

	public ECParameterSpec p256ParameterSpec() {
		ECField javaField = new ECFieldFp(fp256.getNumberOfElements());
		java.security.spec.EllipticCurve curveSpec = new java.security.spec.EllipticCurve(javaField, p256A.getValue(),
				p256B.getValue());
		ECPoint g = toECPoint(p256Generator);
		ECParameterSpec ecSpec = new ECParameterSpec(curveSpec, g, p256Order, 1);
		return ecSpec;
	}

	public ProjectivePoint<PFE> fromECPoint(ECPoint point) {
		if (point.equals(ECPoint.POINT_INFINITY)) {
			return p256.neutral();
		}
		return new ProjectivePoint<PFE>(fp256, fp256.getElement(point.getAffineX()),
				fp256.getElement(point.getAffineY()), fp256.one());
	}

	public ECPoint toECPoint(ProjectivePoint<PFE> point) {
		if (point.equals(p256.neutral())) {
			return ECPoint.POINT_INFINITY;
		}
		return new ECPoint(point.getDehomogenisedCoord(1, 3).getValue(), point.getDehomogenisedCoord(2, 3).getValue());
	}

	public KeyPair makeECKeysJava() {
		try {
			KeyPairGenerator keyGen = KeyPairGenerator.getInstance("EC");
			keyGen.initialize(p256ParameterSpec());
			return keyGen.generateKeyPair();
		} catch (Exception e) {
			return null;
		}
	}

	public byte[] computeSharedSecret(BigInteger privateKey, ProjectivePoint<PFE> publicKey) {
		ProjectivePoint<PFE> sharedSecretPoint = p256.multiply(privateKey, publicKey);
		System.out.println(sharedSecretPoint);
		byte[] sharedSecret = sharedSecretPoint.getDehomogenisedCoord(1, 3).getValue().toByteArray();
		if (sharedSecret[0] == 0) {
			return Arrays.copyOfRange(sharedSecret, 1, sharedSecret.length);
		} else {
			return sharedSecret;
		}
	}

	public byte[] computeSharedSecretJava(BigInteger privateKey, ECPoint publicKey) {
		try {
			KeyFactory kf = KeyFactory.getInstance("EC");
			ECPrivateKeySpec privateKeySpec = new ECPrivateKeySpec(privateKey, p256ParameterSpec());
			ECPublicKeySpec publicKeySpec = new ECPublicKeySpec(publicKey, p256ParameterSpec());

			PrivateKey javaPrivateKey = kf.generatePrivate(privateKeySpec);
			PublicKey javaPublicKey = kf.generatePublic(publicKeySpec);
			KeyAgreement ka = KeyAgreement.getInstance("ECDH");
			ka.init(javaPrivateKey);
			ka.doPhase(javaPublicKey, true /* lastPhase */);
			byte[] secret = ka.generateSecret();
			return secret;
		} catch (Exception ex) {
			System.err.println(ex.toString());
			return new byte[0];
		}
	}
}
