package cryptography;

import cryptography.interfaces.DhScheme;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.interfaces.Element;
import fields.interfaces.Field;
import varieties.curves.elliptic.EllipticCurve;
import varieties.projective.ProjectivePoint;

public class Ecdh<T extends Element<T>> implements DhScheme<ProjectivePoint<T>, ProjectivePoint<T>, IntE, Ecdh<T>> {
	private Field<T> field;
	private EllipticCurve<T> curve;
	private ProjectivePoint<T> generator;
	private IntE order;

	public Ecdh(Field<T> field, EllipticCurve<T> curve, ProjectivePoint<T> generator) {
		this.field = field;
		this.curve = curve;
		this.generator = generator;
		this.order = Integers.z().getInteger(curve.getOrder(generator));
	}

	@Override
	public ProjectivePoint<T> createPublicKey(IntE privateKey) {
		return curve.multiply(privateKey, generator);
	}

	@Override
	public IntE createPrivateKey(Role role) {
		return Integers.z().getRandomElement(order);
	}

	@Override
	public ProjectivePoint<T> agree(IntE privateKey, ProjectivePoint<T> publicKey) {
		if (!curve.hasRationalPoint(publicKey)) {
			throw new ArithmeticException("Point not on curve!");
		}
		if (!curve.weilPairing(publicKey, generator).equals(field.one())) {
			throw new ArithmeticException("Point not in span!");
		}
		return curve.multiply(privateKey, publicKey);
	}

}
