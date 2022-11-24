package cryptography;

import java.math.BigInteger;

import cryptography.SchnorrSignature.Signature;
import cryptography.interfaces.Encoder;
import cryptography.interfaces.HashFunction;
import cryptography.interfaces.SignatureScheme;
import fields.finitefields.PrimeField;
import fields.finitefields.PrimeField.PFE;
import fields.helper.AbstractElement;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.interfaces.Element;
import fields.interfaces.Group;
import varieties.curves.elliptic.EllipticCurve;
import varieties.projective.ProjectivePoint;

public class SchnorrSignature<T extends Element<T>> implements SignatureScheme<Signature, T, PFE, SchnorrSignature<T>> {
	public static class Signature extends AbstractElement<Signature> {
		private PFE e;
		private PFE s;

		private Signature(PFE s, PFE e) {
			this.e = e;
			this.s = s;
		}

		@Override
		public int compareTo(Signature o) {
			int cmp = e.compareTo(o.e);
			if (cmp != 0) {
				return cmp;
			}
			return s.compareTo(o.s);
		}
	}

	public static <T extends Element<T>> SchnorrSignature<ProjectivePoint<T>> ecSchnorr(HashFunction hash,
			Encoder<T> encoder, EllipticCurve<T> curve, ProjectivePoint<T> generator) {
		return new SchnorrSignature<>(
				Encoder.fromFunction((ProjectivePoint<T> t) -> encoder.encode(t.getDehomogenisedCoord(1, 3))), hash,
				curve, generator);
	}

	private Encoder<T> encoder;
	private HashFunction hash;
	private PrimeField field;
	private Group<T> group;
	private T generator;
	private IntE order;

	public SchnorrSignature(Encoder<T> encoder, HashFunction hash, Group<T> group, T generator) {
		this.encoder = encoder;
		this.hash = hash;
		this.group = group;
		this.generator = generator;
		this.order = Integers.z().getInteger(group.getOrder(generator));
		this.field = PrimeField.getPrimeField(order);
	}

	@Override
	public T createPublicKey(PFE privateKey) {
		return group.power(privateKey.getValue(), generator);
	}

	@Override
	public PFE createPrivateKey() {
		return field.getRandomElement();
	}

	@Override
	public Signature sign(byte[] message, PFE privateKey) {
		PFE nonce = field.getRandomElement();
		T r = group.power(nonce.getValue(), generator);
		byte[] encoded = encoder.encode(r);
		PFE e = field.getInteger(new BigInteger(hash.evaluate(ByteArray.concatenate(encoded, message).array())));
		PFE s = field.subtract(nonce, field.multiply(e, privateKey));
		return new Signature(s, e);
	}

	@Override
	public boolean verify(byte[] message, Signature signature, T publicKey) {
		T reconstructed = group.operate(group.power(signature.s.getValue(), generator),
				group.power(signature.e.getValue(), publicKey));
		byte[] encoded = encoder.encode(reconstructed);
		PFE e = field.getInteger(new BigInteger(hash.evaluate(ByteArray.concatenate(encoded, message).array())));
		return e.equals(signature.e);
	}

}
