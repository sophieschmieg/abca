package main;

import java.math.BigInteger;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import fields.interfaces.Element;
import fields.interfaces.Field;
import util.MiscAlgorithms;
import varieties.ProjectivePoint;
import varieties.curves.EllipticCurve;

public class TestFindRelatedCurves<T extends Element<T>> {
	private EllipticCurve<T> curve;
	private ProjectivePoint<T> generator;
	private List<RelatedCurveAttack<T>> rcas = new ArrayList<>();

	public TestFindRelatedCurves(Field<T> field, T a, T b, ProjectivePoint<T> generator) {
		this.curve = new EllipticCurve<>(field, a, b);
		this.generator = generator;
		BigInteger q = field.getNumberOfElements();
		BigInteger number = curve.getNumberOfElements();
		BigInteger target = BigInteger.ONE;
		BigInteger max = BigInteger.valueOf(q.bitLength()).add(BigInteger.TEN);
		System.out.println(max);
		Map<BigInteger, RelatedCurveAttack<T>> rcas = new TreeMap<>();
		for (T bVariant : field) {
			if (target.compareTo(number) > 0) {
				break;
			}
			EllipticCurve<T> curve;
			try {
				curve = new EllipticCurve<>(field, a, bVariant);
			} catch (ArithmeticException e) {
				continue;
			}
			System.out.println("Testing " + curve);
			BigInteger numPoints = curve.getNumberOfElements();
			System.out.println("Number of points: " + numPoints);

			Map<BigInteger, Integer> decomp = MiscAlgorithms.naivePrimeDecomposition(numPoints, max);
			for (BigInteger prime : decomp.keySet()) {
				BigInteger primePower = prime;
				int power = 1;
				while (power <= decomp.get(prime) && primePower.multiply(prime).compareTo(max) <= 0) {
					power++;
					primePower = primePower.multiply(prime);
				}
				if (rcas.containsKey(prime) && power <= rcas.get(prime).power) {
					continue;
				}
				System.out.println("Checking for " + prime + "^" + power + " torsion points.");
				List<ProjectivePoint<T>> torsionPoints = curve.getTorsionPoints(primePower.intValueExact());
				int maxPower = 0;
				ProjectivePoint<T> torsionPoint = curve.neutral();
				for (ProjectivePoint<T> point : torsionPoints) {
					for (int pow = maxPower; pow <= power; pow++) {
						if (curve.multiply(prime.pow(pow), point).equals(curve.neutral())) {
							if (pow > maxPower) {
								maxPower = pow;
								torsionPoint = point;
							}
							break;
						}
					}
				}
				System.out.println("Found non trivial " + prime + "^" + maxPower + " torsion point on " + curve);
				if (!rcas.containsKey(prime) || rcas.get(prime).power < maxPower) {
					if (rcas.containsKey(prime)) {
						RelatedCurveAttack<T> rca = rcas.get(prime);
						target.divide(rca.prime.pow(rca.power));
					}
					rcas.put(prime, new RelatedCurveAttack<>(bVariant, prime, maxPower, torsionPoint, curve));
					target = target.multiply(prime.pow(maxPower));
					System.out.println("New target is " + target);
				}
			}
		}
		this.rcas = new ArrayList<>();
		this.rcas.addAll(rcas.values());
		System.out.println("Attack vectors: " + this.rcas);
	}

	public void doAttack() {
		Oracle<T> oracle = new Oracle<>(curve, generator);
		List<BigInteger> moduli = new ArrayList<>();
		List<BigInteger> moduloValues = new ArrayList<>();
		for (RelatedCurveAttack<T> rca : rcas) {
			BigInteger prime = rca.prime;
			int power = rca.power;
			BigInteger primePower = prime.pow(power);
			moduli.add(primePower);
			for (ProjectivePoint<T> point : rca.torsionPoints.keySet()) {
				if (oracle.oracle(rca.torsionPoint, point)) {
					moduloValues.add(rca.torsionPoints.get(point));
					break;
				}
			}
		}
		System.out.println(moduli);
		System.out.println(moduloValues);
		BigInteger guess = MiscAlgorithms.chineseRemainder(moduloValues, moduli, false);
		guess = guess.mod(this.curve.getNumberOfElements());
		System.out.println("Guess for private key: " + guess);
		System.out.println("Number of oracle invocations: " + oracle.getGuessCounter());
		System.out.println("Guess correct: " + oracle.guessPrivateKey(guess));
	}

	private static class RelatedCurveAttack<T extends Element<T>> {
		private T b;
		private BigInteger prime;
		private int power;
		private ProjectivePoint<T> torsionPoint;
		private Map<ProjectivePoint<T>, BigInteger> torsionPoints;

		public RelatedCurveAttack(T b, BigInteger prime, int power, ProjectivePoint<T> torsionPoint,
				EllipticCurve<T> curve) {
			this.b = b;
			this.prime = prime;
			this.power = power;
			this.torsionPoint = torsionPoint;
			this.torsionPoints = new TreeMap<>();
			for (BigInteger i = BigInteger.ZERO; i.compareTo(prime.pow(power)) < 0; i = i.add(BigInteger.ONE)) {
				this.torsionPoints.put(curve.multiply(i, torsionPoint), i);
			}
		}

		@Override
		public String toString() {
			return "(b: " + this.b + ", [" + this.prime + "^" + this.power + "]" + this.torsionPoint + ")";
		}
	}

	public static class Oracle<T extends Element<T>> {
		private BigInteger privateKey;
		private ProjectivePoint<T> publicKey;
		private EllipticCurve<T> curve;
		private int guessCounter;

		public Oracle(EllipticCurve<T> curve, ProjectivePoint<T> generator) {
			this.curve = curve;
			this.guessCounter = 0;
			BigInteger order = curve.getNumberOfElements();
			SecureRandom rand = new SecureRandom();
			byte[] privateBytes = new byte[order.toByteArray().length];
			do {
				rand.nextBytes(privateBytes);
				this.privateKey = new BigInteger(privateBytes);
			} while (privateKey.compareTo(order) >= 0 || privateKey.compareTo(BigInteger.ZERO) < 0);
			this.publicKey = curve.multiply(privateKey, generator);
		}

		public ProjectivePoint<T> getPublicKey() {
			return publicKey;
		}

		public boolean oracle(ProjectivePoint<T> point, ProjectivePoint<T> sharedSecret) {
			guessCounter++;
			return curve.multiply(privateKey, point).equals(sharedSecret);
		}

		public int getGuessCounter() {
			return guessCounter;
		}

		public boolean guessPrivateKey(BigInteger guess) {
			return guess.equals(privateKey);
		}
	}
}
