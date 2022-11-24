package cryptography;

import java.util.function.Function;

import cryptography.interfaces.DhScheme;
import cryptography.interfaces.Encoder;
import cryptography.interfaces.ExtendedOutputFunction;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.interfaces.Element;
import fields.interfaces.Group;
import varieties.curves.elliptic.EllipticCurve;
import varieties.projective.ProjectivePoint;

public class DiffieHellman<T extends Element<T>> implements DhScheme<T, IntE, DiffieHellman<T>> {
	public static <T extends Element<T>> DiffieHellman<ProjectivePoint<T>> ecdh(ExtendedOutputFunction kdf,
			Encoder<T> encoder, EllipticCurve<T> curve, ProjectivePoint<T> generator) {
		return new DiffieHellman<>(kdf,
				Encoder.fromFunction((ProjectivePoint<T> t) -> encoder.encode(t.getDehomogenisedCoord(1, 3))),
				(ProjectivePoint<T> publicKey) -> {
					if (!curve.hasRationalPoint(publicKey)) {
						return false;
					}
					if (!curve.weilPairing(publicKey, generator).equals(curve.getField().one())) {
						return false;
					}
					return true;
				}, curve, generator);

	}

	private ExtendedOutputFunction kdf;
	private Encoder<T> encoder;
	private Function<T, Boolean> publicKeyCheck;
	private Group<T> group;
	private T generator;
	private IntE order;

	public DiffieHellman(ExtendedOutputFunction kdf, Encoder<T> encoder, Group<T> group, T generator) {
		this(kdf, encoder, (T t) -> true, group, generator);
	}

	public DiffieHellman(ExtendedOutputFunction kdf, Encoder<T> encoder, Function<T, Boolean> publicKeyCheck,
			Group<T> group, T generator) {
		this.kdf = kdf;
		this.encoder = encoder;
		this.publicKeyCheck = publicKeyCheck;
		this.group = group;
		this.generator = generator;
		this.order = Integers.z().getInteger(group.getOrder(generator));
	}

	@Override
	public T createPublicKey(IntE privateKey) {
		return group.power(privateKey.getValue(), generator);
	}

	@Override
	public IntE createPrivateKey(Role role) {
		return Integers.z().getRandomElement(order);
	}

	@Override
	public VariableLengthKey agree(IntE privateKey, T publicKey) {
		if (!publicKeyCheck.apply(publicKey)) {
			throw new ArithmeticException("Invalid public key!");
		}
		return new VariableLengthKey(encoder.encode(group.power(privateKey.getValue(), publicKey)), kdf);
	}

}
