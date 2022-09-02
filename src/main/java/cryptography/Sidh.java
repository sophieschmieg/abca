package cryptography;

import java.math.BigInteger;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.List;

import cryptography.Sidh.SidhPrivateKey;
import cryptography.Sidh.SidhPublicKey;
import cryptography.interfaces.DhScheme;
import fields.finitefields.FiniteField;
import fields.finitefields.FiniteField.FFE;
import fields.helper.AbstractElement;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import util.MiscAlgorithms;
import varieties.curves.elliptic.EllipticCurve;
import varieties.curves.elliptic.Isogeny;
import varieties.curves.elliptic.KernelIsogeny;
import varieties.projective.ProjectivePoint;

public class Sidh implements DhScheme<FFE, SidhPublicKey, SidhPrivateKey, Sidh> {
	private FiniteField field;
	private EllipticCurve<FFE> publicCurve;
	private List<ProjectivePoint<FFE>> aliceBasis;
	private BigInteger aliceOrder;
	private List<ProjectivePoint<FFE>> bobBasis;
	private BigInteger bobOrder;

	public Sidh(IntE prime) {
		Integers z = Integers.z();
		if (!z.remainder(prime, z.getInteger(4)).equals(z.getInteger(3)) || !z.isPrime(prime)) {
			throw new ArithmeticException("Need a prime where -1 is not a square!");
		}
		this.field = FiniteField.getFiniteField(prime.getValue(), 2);
		this.publicCurve = EllipticCurve.fromJInvariant(field, field.getInteger(1728));
		BigInteger one = BigInteger.ONE;
		BigInteger two = BigInteger.TWO;
		BigInteger three = BigInteger.valueOf(3);
		BigInteger order = publicCurve.getField().characteristic().add(one);
		this.aliceOrder = getPrimeTorsionOrder(two, order);
		this.bobOrder = getPrimeTorsionOrder(three, order);
		this.aliceBasis = publicCurve.getTorsionPointBasis(aliceOrder);
		this.bobBasis = publicCurve.getTorsionPointBasis(bobOrder);
	}

	private BigInteger getPrimeTorsionOrder(BigInteger prime, BigInteger order) {
		BigInteger zero = BigInteger.ZERO;
		BigInteger one = BigInteger.ONE;
		BigInteger primeTorsionPointOrder = one;
		while (order.mod(primeTorsionPointOrder.multiply(prime)).equals(zero)) {
			primeTorsionPointOrder = primeTorsionPointOrder.multiply(prime);
		}
		if (primeTorsionPointOrder.equals(one)) {
			throw new ArithmeticException("Not enough " + prime + " torsion points");
		}
		return primeTorsionPointOrder;

	}

	@Override
	public SidhPrivateKey createPrivateKey(Role role) {
		BigInteger order = role.equals(Role.Alice) ? aliceOrder : bobOrder;
		List<ProjectivePoint<FFE>> basis = role.equals(Role.Alice) ? aliceBasis : bobBasis;
		List<ProjectivePoint<FFE>> partnerBasis = role.equals(Role.Alice) ? bobBasis : aliceBasis;
		List<BigInteger> key = new ArrayList<BigInteger>();
		SecureRandom rand = new SecureRandom();
		boolean usePoint1 = rand.nextBoolean();
		BigInteger privateKey = MiscAlgorithms.randomBigInteger(rand, order);
		key.add(usePoint1 ? privateKey : BigInteger.ONE);
		key.add(usePoint1 ? BigInteger.ONE : privateKey);
		return new SidhPrivateKey(key, order, basis, partnerBasis);
	}
	
	@Override
	public SidhPublicKey createPublicKey(SidhPrivateKey privateKey) {
		return privateKey.getPublicKey();
	}
	
	@Override
	public FFE agree(SidhPrivateKey privateKey, SidhPublicKey publicKey) {
		return privateKey.computeSharedSecret(publicKey);
	}
	
	public FiniteField getField() {
		return field;
	}

	public EllipticCurve<FFE> getPublicCurve() {
		return publicCurve;
	}

	public List<ProjectivePoint<FFE>> getAliceBasis() {
		return aliceBasis;
	}

	public BigInteger getAliceOrder() {
		return aliceOrder;
	}

	public List<ProjectivePoint<FFE>> getBobBasis() {
		return bobBasis;
	}

	public BigInteger getBobOrder() {
		return bobOrder;
	}

	public String toString() {
		return publicCurve.toString() + " / " + publicCurve.getField().toString() + ", (" + aliceBasis.get(0) + ", "
				+ aliceBasis.get(1) + "), (" + bobBasis.get(0) + ", " + bobBasis.get(1) + ")";
	}

	public class SidhPrivateKey extends AbstractElement<SidhPrivateKey> {
		private List<ProjectivePoint<FFE>> basis;
		private BigInteger order;
		private List<ProjectivePoint<FFE>> partnerBasis;
		private List<BigInteger> coefficients;
		private Isogeny<FFE> isogeny;

		public SidhPrivateKey(List<BigInteger> coefficients, BigInteger order, List<ProjectivePoint<FFE>> basis,
				List<ProjectivePoint<FFE>> partnerBasis) {
			this.coefficients = coefficients;
			ProjectivePoint<FFE> point1 = publicCurve.multiply(coefficients.get(0), basis.get(0));
			ProjectivePoint<FFE> point2 = publicCurve.multiply(coefficients.get(1), basis.get(1));
			ProjectivePoint<FFE> point = publicCurve.add(point1, point2);
			this.isogeny = new KernelIsogeny<>(publicCurve, point, order);
			this.basis = basis;
			this.order = order;
			this.partnerBasis = partnerBasis;
		}

		private SidhPublicKey getPublicKey() {
			List<ProjectivePoint<FFE>> inducedBasis = new ArrayList<>();
			for (ProjectivePoint<FFE> point : partnerBasis) {
				inducedBasis.add(isogeny.evaluate(point));
			}
			return new SidhPublicKey(isogeny.getRange(), inducedBasis);
		}

		private FFE computeSharedSecret(SidhPublicKey partnerKey) {
			ProjectivePoint<FFE> point1 = partnerKey.curve.multiply(coefficients.get(0),
					partnerKey.inducedBasis.get(0));
			ProjectivePoint<FFE> point2 = partnerKey.curve.multiply(coefficients.get(1),
					partnerKey.inducedBasis.get(1));
			ProjectivePoint<FFE> point = partnerKey.curve.add(point1, point2);
			Isogeny<FFE> isogeny = new KernelIsogeny<>(partnerKey.curve, point, order);
			return isogeny.getRange().jInvariant();
		}

		public String toString() {
			return "[" + coefficients.get(0) + "]" + basis.get(0) + " + " + "[" + coefficients.get(1) + "]"
					+ basis.get(1);
		}

		@Override
		public int compareTo(SidhPrivateKey o) {
			return isogeny.compareTo(o.isogeny);
		}
	}

	public class SidhPublicKey extends AbstractElement<SidhPublicKey> {
		private EllipticCurve<FFE> curve;
		private List<ProjectivePoint<FFE>> inducedBasis;

		public SidhPublicKey(EllipticCurve<FFE> curve, List<ProjectivePoint<FFE>> inducedBasis) {
			this.curve = curve;
			this.inducedBasis = inducedBasis;
		}

		public String toString() {
			return this.curve + " (" + inducedBasis.get(0) + ", " + inducedBasis.get(1) + ")";
		}

		@Override
		public int compareTo(SidhPublicKey o) {
			return curve.jInvariant().compareTo(o.curve.jInvariant());
		}
	}
}
