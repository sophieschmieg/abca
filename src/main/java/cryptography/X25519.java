package cryptography;

import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.interfaces.MathMap;
import varieties.curves.elliptic.EllipticCurve;
import varieties.projective.ProjectivePoint;

public class X25519 implements MathMap<ByteArray, ByteArray> {
	private static PrimeField getField() {
		Integers z = Integers.z();
		IntE prime = z.subtract(z.power(z.getInteger(2), 255), z.getInteger(19));
		return PrimeField.getPrimeField(prime);
	}

	private final static PrimeField field = getField();
	private final static EllipticCurve<PFE> curve25519 = new EllipticCurve<>(field, field.zero(),
			field.getInteger(486662), field.zero(), field.one(), field.zero());
	private final static EllipticCurve<PFE> twist = curve25519.getQuadraticTwist();
	private final static PFE twistMultiplier = field.divide(twist.getA2(), curve25519.getA2());
	public final static ByteArray GENERATOR = new ByteArray(new byte[] { 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 });
	private IntE key;

	public X25519(byte[] key) {
		if (key.length != 32) {
			throw new ArithmeticException("key invalid");
		}
		key[0] &= 248;
		key[31] &= 127;
		key[31] |= 64;
		this.key = decodeLittleEndian(key);
	}

	private IntE decodeLittleEndian(byte[] t) {
		Integers z = Integers.z();
		IntE result = z.zero();
		for (int i = t.length - 1; i >= 0; i--) {
			result = z.multiply(256, result);
			result = z.add(z.getInteger(((int) t[i]) & 255), result);
		}
		return result;
	}

	@Override
	public ByteArray evaluate(ByteArray t) {
		byte[] array = t.array();
		array[31] &= 127;
		PFE u = field.getInteger(decodeLittleEndian(array));
		PFE v2 = field.add(field.power(u, 3), field.multiply(curve25519.getA2(), u, u), u);
		PFE shared;
		if (field.hasSqrt(v2)) {
			PFE v = field.sqrt(v2).keySet().iterator().next();
			ProjectivePoint<PFE> pk = new ProjectivePoint<>(field, u, v, field.one());
			shared = curve25519.multiply(key, pk).getDehomogenisedCoord(1, 3);
		} else {
			u = field.multiply(twistMultiplier, u);
			v2 = field.multiply(field.power(twistMultiplier, 3), v2);
			PFE v = field.sqrt(v2).keySet().iterator().next();
			ProjectivePoint<PFE> pk = new ProjectivePoint<>(field, u, v, field.one());
			shared = field.divide(twist.multiply(key, pk).getDehomogenisedCoord(1, 3), twistMultiplier);
		}
		Integers z = Integers.z();
		IntE twoFiveSix = z.getInteger(256);
		IntE x = z.lift(shared);
		byte[] enc = new byte[32];
		for (int i = 0; i < enc.length; i++) {
			enc[i] = (byte) z.remainder(x, twoFiveSix).intValueExact();
			x = z.divide(x, twoFiveSix);
		}
		return new ByteArray(enc);
	}
}
