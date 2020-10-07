package main;

import java.math.BigInteger;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.List;

import fields.interfaces.Element;
import fields.interfaces.Field;
import util.MiscAlgorithms;
import varieties.ProjectivePoint;
import varieties.curves.EllipticCurve;
import varieties.curves.Isogeny;
import varieties.curves.KernelIsogeny;

public class TestSike<T extends Element<T>> {
	private Field<T> field;
	private EllipticCurve<T> supersingular;

	public TestSike(Field<T> field) {
		if (!field.characteristic().multiply(field.characteristic()).equals(field.getNumberOfElements())) {
			throw new ArithmeticException("Field needs to be a quadratic extension!");
		}
		this.field = field;
		this.supersingular = new EllipticCurve<>(field, field.negative(field.one()), field.zero());
		System.out.println(supersingular + "; j = " + supersingular.jInvariant());
		System.out.println("Number of points: " + supersingular.getNumberOfElements());
		SidhParameters parameters = makeSidhParameters();
		System.out.println("Parameters: " + parameters);
		SidhPrivateKey alicePrivateKey = makePrivateKey(parameters);
		System.out.println("Private key Alice: " + alicePrivateKey);
		SidhPublicKey alicePublicKey = alicePrivateKey.getPublicKey();
		System.out.println("Public key Alice: " + alicePublicKey);
		SidhPrivateKey bobPrivateKey = makePrivateKey(parameters.swap());
		System.out.println("Private key Bob: " + bobPrivateKey);
		SidhPublicKey bobPublicKey = bobPrivateKey.getPublicKey();
		System.out.println("Public key Bob: " + bobPublicKey);

		System.out.println("shared secret Alice: " + alicePrivateKey.computeSharedSecret(bobPublicKey));
		System.out.println("shared secret Bob: " + bobPrivateKey.computeSharedSecret(alicePublicKey));
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

	private List<ProjectivePoint<T>> computeBasis(BigInteger prime, BigInteger order, BigInteger basisOrder) {
		List<ProjectivePoint<T>> basis = new ArrayList<>();
		ProjectivePoint<T> basis1;
		do {
			System.out.println("Count basis1 up");

			basis1 = supersingular.multiply(order.divide(basisOrder), supersingular.getRandomElement());
		} while (!supersingular.multiply(basisOrder, basis1).equals(supersingular.neutral())
				|| supersingular.multiply(basisOrder.divide(prime), basis1).equals(supersingular.neutral()));
		ProjectivePoint<T> basis2;
		do {
			System.out.println("Count basis2 up");
			basis2 = supersingular.multiply(order.divide(basisOrder), supersingular.getRandomElement());
		} while (!supersingular.multiply(basisOrder, basis2).equals(supersingular.neutral())
				|| supersingular.multiply(basisOrder.divide(prime), basis2).equals(supersingular.neutral())
				|| supersingular.weilPairing(basisOrder, basis1, basis2).equals(field.one()));
		basis.add(basis1);
		basis.add(basis2);
		return basis;
	}

	private SidhParameters makeSidhParameters() {
		BigInteger one = BigInteger.ONE;
		BigInteger two = BigInteger.TWO;
		BigInteger three = BigInteger.valueOf(3);
		BigInteger order = supersingular.getField().characteristic().add(one);
		BigInteger twoPowerOrder = getPrimeTorsionOrder(two, order);
		BigInteger threePowerOrder = getPrimeTorsionOrder(three, order);
		List<ProjectivePoint<T>> twoBasis = computeBasis(two, order, twoPowerOrder);
		List<ProjectivePoint<T>> threeBasis = computeBasis(three, order, threePowerOrder);
		return new SidhParameters(supersingular, twoBasis, twoPowerOrder, threeBasis, threePowerOrder);
	}

	private SidhPrivateKey makePrivateKey(SidhParameters parameters) {
		List<BigInteger> key = new ArrayList<BigInteger>();
		SecureRandom rand = new SecureRandom();
		boolean usePoint1 = rand.nextBoolean();
		BigInteger privateKey = MiscAlgorithms.randomBigInteger(rand, parameters.order);
		key.add(usePoint1 ? privateKey : BigInteger.ONE);
		key.add(usePoint1 ? BigInteger.ONE : privateKey);
		return new SidhPrivateKey(parameters, key);
	}

	public class SidhParameters {
		private EllipticCurve<T> publicCurve;
		private List<ProjectivePoint<T>> basis;
		private BigInteger order;
		private List<ProjectivePoint<T>> partnerBasis;
		private BigInteger partnerOrder;

		public SidhParameters(EllipticCurve<T> publicCurve, List<ProjectivePoint<T>> basis, BigInteger order,
				List<ProjectivePoint<T>> partnerBasis, BigInteger partnerOrder) {
			this.publicCurve = publicCurve;
			this.basis = basis;
			this.order = order;
			this.partnerBasis = partnerBasis;
			this.partnerOrder = partnerOrder;
		}

		public SidhParameters swap() {
			return new SidhParameters(publicCurve, partnerBasis, partnerOrder, basis, order);
		}

		@SuppressWarnings("unchecked")
		public boolean equals(Object o) {
			if (!(o instanceof TestSike.SidhParameters)) {
				return false;
			}
			SidhParameters p = (SidhParameters) o;
			return publicCurve.equals(p.publicCurve) && basis.equals(p.basis) && order.equals(p.order)
					&& partnerBasis.equals(p.partnerBasis) && partnerOrder.equals(p.partnerOrder);
		}

		public String toString() {
			return publicCurve.toString() + " / " + publicCurve.getField().toString() + ", (" + basis.get(0) + ", "
					+ basis.get(1) + "), (" + partnerBasis.get(0) + ", " + partnerBasis.get(1) + ")";
		}
	}

	public class SidhPrivateKey {
		private SidhParameters parameters;
		private List<BigInteger> coefficients;
		private Isogeny<T> isogeny;

		public SidhPrivateKey(SidhParameters parameters, List<BigInteger> coefficients) {
			this.parameters = parameters;
			this.coefficients = coefficients;
			ProjectivePoint<T> point1 = parameters.publicCurve.multiply(coefficients.get(0), parameters.basis.get(0));
			ProjectivePoint<T> point2 = parameters.publicCurve.multiply(coefficients.get(1), parameters.basis.get(1));
			ProjectivePoint<T> point = parameters.publicCurve.add(point1, point2);
			this.isogeny = new KernelIsogeny<T>(parameters.publicCurve, point, parameters.order);
		}

		public SidhPublicKey getPublicKey() {
			List<ProjectivePoint<T>> inducedBasis = new ArrayList<>();
			for (ProjectivePoint<T> point : parameters.partnerBasis) {
				inducedBasis.add(isogeny.evaluate(point));
			}
			return new SidhPublicKey(parameters, isogeny.getRange(), inducedBasis);
		}

		public T computeSharedSecret(SidhPublicKey partnerKey) {
			if (!partnerKey.parameters.swap().equals(this.parameters)) {
				throw new RuntimeException("Parameters do not match!");
			}
			ProjectivePoint<T> point1 = partnerKey.curve.multiply(coefficients.get(0), partnerKey.inducedBasis.get(0));
			ProjectivePoint<T> point2 = partnerKey.curve.multiply(coefficients.get(1), partnerKey.inducedBasis.get(1));
			ProjectivePoint<T> point = partnerKey.curve.add(point1, point2);
			Isogeny<T> isogeny = new KernelIsogeny<T>(partnerKey.curve, point, parameters.order);
			return isogeny.getRange().jInvariant();
		}

		public String toString() {
			return "[" + coefficients.get(0) + "]" + parameters.basis.get(0) + " + " + "[" + coefficients.get(1) + "]"
					+ parameters.basis.get(1);
		}
	}

	public class SidhPublicKey {
		private SidhParameters parameters;
		private EllipticCurve<T> curve;
		private List<ProjectivePoint<T>> inducedBasis;

		public SidhPublicKey(SidhParameters parameters, EllipticCurve<T> curve, List<ProjectivePoint<T>> inducedBasis) {
			this.parameters = parameters;
			this.curve = curve;
			this.inducedBasis = inducedBasis;
		}

		public String toString() {
			return this.curve + " (" + inducedBasis.get(0) + ", " + inducedBasis.get(1) + ")";
		}
	}
}
