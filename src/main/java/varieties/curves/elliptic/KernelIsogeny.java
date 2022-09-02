package varieties.curves.elliptic;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Deque;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import fields.helper.AbstractElement;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.interfaces.Element;
import fields.interfaces.Ring.FactorizationResult;
import varieties.projective.ProjectiveMorphism;
import varieties.projective.ProjectivePoint;

public class KernelIsogeny<T extends Element<T>> extends AbstractElement<Isogeny<T>> implements Isogeny<T> {
	private Deque<Isogeny<T>> isogenies;
	private IntE degree;
	private List<ProjectivePoint<T>> kernelGenerators;

	public KernelIsogeny(EllipticCurve<T> domain, ProjectivePoint<T> kernel) {
		this(domain, Collections.singletonMap(kernel, domain.getNumberOfElements()));
	}

	public KernelIsogeny(EllipticCurve<T> domain, ProjectivePoint<T> kernel, BigInteger order) {
		this(domain, Collections.singletonMap(kernel, order));
	}

	private static <T extends Element<T>> Map<ProjectivePoint<T>, BigInteger> fillMap(BigInteger order,
			Set<ProjectivePoint<T>> kernel) {
		Map<ProjectivePoint<T>, BigInteger> map = new TreeMap<>();
		for (ProjectivePoint<T> p : kernel) {
			map.put(p, order);
		}
		return map;
	}

	public KernelIsogeny(EllipticCurve<T> domain, Set<ProjectivePoint<T>> kernel) {
		this(domain, fillMap(domain.getNumberOfElements(), kernel));
	}

	public KernelIsogeny(EllipticCurve<T> domain, Set<ProjectivePoint<T>> kernel, BigInteger order) {
		this(domain, fillMap(order, kernel));
	}

	public KernelIsogeny(EllipticCurve<T> domain, Map<ProjectivePoint<T>, BigInteger> kernel) {
		Integers z = Integers.z();
		this.isogenies = new LinkedList<>();
		this.degree = z.one();
		Map<ProjectivePoint<T>, IntE> simplified = simplifyTorsionPoints(domain, kernel);
		this.kernelGenerators = new ArrayList<>();
		for (ProjectivePoint<T> point : simplified.keySet()) {
			kernelGenerators.add(point);
			degree = z.multiply(isogenyLadder(domain, point, simplified.get(point), isogenies), degree);
		}
//		EllipticCurve<T> currentRange = domain;
//		for (ProjectivePoint<T> kernelGenerator : kernel.keySet()) {
//			BigInteger order = kernel.get(kernelGenerator);
//			if (!domain.hasRationalPoint(kernelGenerator)) {
//				throw new ArithmeticException("Point " + kernelGenerator + " not on curve!");
//			}
//			if (!domain.multiply(order, kernelGenerator).equals(domain.neutral())) {
//				throw new ArithmeticException("Order hint " + order + " for " + kernelGenerator + " was wrong!");
//			}
//			BigInteger prime = BigInteger.TWO;
//			Map<BigInteger, Integer> orderPrimeDecomposition = new TreeMap<>();
//			BigInteger reducedOrder = order;
//			while (!reducedOrder.equals(BigInteger.ONE)) {
//				if (!reducedOrder.mod(prime).equals(BigInteger.ZERO)) {
//					prime = prime.nextProbablePrime();
//					continue;
//				}
//				reducedOrder = order.divide(prime);
//				if (domain.multiply(reducedOrder, kernelGenerator).equals(domain.neutral())) {
//					order = reducedOrder;
//					continue;
//				}
//				int i = 1;
//				while (reducedOrder.mod(prime).equals(BigInteger.ZERO)) {
//					i++;
//					reducedOrder = reducedOrder.divide(prime);
//				}
//				orderPrimeDecomposition.put(prime, i);
//				prime = prime.nextProbablePrime();
//			}
//			ProjectivePoint<T> kernelGeneratorImage = this.evaluate(kernelGenerator);
//			for (BigInteger orderPrime : orderPrimeDecomposition.keySet()) {
//				int exponent = orderPrimeDecomposition.get(orderPrime);
//				BigInteger primePower = orderPrime.pow(exponent);
//				while (!primePower.equals(BigInteger.ONE)) {
//					primePower = primePower.divide(orderPrime);
//					ProjectivePoint<T> powerKernelPoint = currentRange.multiply(primePower, kernelGeneratorImage);
//					Isogeny<T> nextIsogeny = new KernelPointIsogeny<T>(currentRange, powerKernelPoint,
//							orderPrime.intValueExact());
//					this.isogenies.add(nextIsogeny);
//					this.degree = this.degree.multiply(nextIsogeny.getDegree());
//					currentRange = nextIsogeny.getRange();
//					kernelGeneratorImage = nextIsogeny.evaluate(kernelGeneratorImage);
//				}
//			}
//		}
	}

	private Map<ProjectivePoint<T>, IntE> simplifyTorsionPoints(EllipticCurve<T> domain,
			Map<ProjectivePoint<T>, BigInteger> points) {
		Integers z = Integers.z();
		Map<IntE, List<ProjectivePoint<T>>> basisPerPrime = new TreeMap<>();
		Map<IntE, List<Integer>> basisOrderExponent = new TreeMap<>();
		for (ProjectivePoint<T> point : points.keySet()) {
			IntE order = z.getInteger(points.get(point));
			FactorizationResult<IntE, IntE> factors = z.uniqueFactorization(order);
			for (IntE prime : factors.primeFactors()) {
				int power = factors.multiplicity(prime);
				IntE multiplier = z.divideChecked(order, z.power(prime, power));
				ProjectivePoint<T> torsion = domain.multiply(multiplier, point);
				while (power > 0 && domain.multiply(z.power(prime, power - 1), torsion).equals(domain.neutral())) {
					power--;
				}
				if (power == 0) {
					continue;
				}
				if (!basisPerPrime.containsKey(prime)) {
					List<ProjectivePoint<T>> basis = new ArrayList<>();
					basis.add(torsion);
					basisPerPrime.put(prime, basis);
					List<Integer> torsionOrder = new ArrayList<>();
					torsionOrder.add(power);
					basisOrderExponent.put(prime, torsionOrder);
					continue;
				}
				List<ProjectivePoint<T>> basis = basisPerPrime.get(prime);
				List<Integer> torsionOrder = basisOrderExponent.get(prime);
				ProjectivePoint<T> highest = power > torsionOrder.get(0) ? torsion : basis.get(0);
				int highestPower = power > torsionOrder.get(0) ? power : torsionOrder.get(0);
				ProjectivePoint<T> second = power > torsionOrder.get(0) ? basis.get(0) : torsion;
				basis.set(0, highest);
				torsionOrder.set(0, highestPower);
				T pairing = domain.weilPairing(highest, second);
				int secondPower = 0;
				while (!pairing.equals(domain.getField().one())) {
					secondPower++;
					pairing = domain.getField().power(pairing, prime);
				}
				if (secondPower == 0) {
					continue;
				}
				if (basis.size() == 1) {
					basis.add(second);
					torsionOrder.add(secondPower);
					continue;
				}
				if (secondPower > torsionOrder.get(1)) {
					basis.set(1, second);
					torsionOrder.set(1, secondPower);
				}
			}
		}
		ProjectivePoint<T> first = domain.neutral();
		ProjectivePoint<T> second = domain.neutral();
		IntE order = z.one();
		for (IntE prime : basisPerPrime.keySet()) {
			List<ProjectivePoint<T>> basis = basisPerPrime.get(prime);
			first = domain.add(first, basis.get(0));
			order = z.multiply(z.power(prime, basisOrderExponent.get(prime).get(0)), order);
			if (basis.size() > 1) {
				second = domain.add(second, basis.get(1));
			}
		}
		Map<ProjectivePoint<T>, IntE> result = new TreeMap<>();
		if (!first.equals(domain.neutral())) {
			result.put(first, order);
		}
		if (!second.equals(domain.neutral())) {
			result.put(second, order);
		}
		return result;
	}

	private IntE isogenyLadder(EllipticCurve<T> domain, ProjectivePoint<T> point, IntE order,
			Deque<Isogeny<T>> ladder) {
		Integers z = Integers.z();
		point = evaluateLadder(point, ladder);
		if (!ladder.isEmpty()) {
			domain = ladder.getLast().getRange();
		}
		FactorizationResult<IntE, IntE> factorization = z.uniqueFactorization(order);
		IntE actualOrder = z.one();
		for (IntE prime : factorization.primeFactors()) {
			IntE multiplier = z.divideChecked(order, z.power(prime, factorization.multiplicity(prime)));
			ProjectivePoint<T> torsion = domain.multiply(multiplier, point);
			int power = factorization.multiplicity(prime);
			while (power > 0 && domain.multiply(z.power(prime, power - 1), torsion).equals(domain.neutral())) {
				power--;
			}
			while (power > 0) {
				ProjectivePoint<T> kernelPoint = domain.multiply(z.power(prime, power - 1), torsion);
				power--;
				actualOrder = z.multiply(prime, actualOrder);
				KernelPointIsogeny<T> isogeny = new KernelPointIsogeny<>(domain, kernelPoint, prime.intValueExact());
				ladder.add(isogeny);
				point = isogeny.evaluate(point);
				torsion = isogeny.evaluate(torsion);
				domain = isogeny.getRange();
			}
		}
		return actualOrder;
	}

	private ProjectivePoint<T> evaluateLadder(ProjectivePoint<T> point, Deque<Isogeny<T>> ladder) {
		for (Isogeny<T> isogeny : ladder) {
			point = isogeny.evaluate(point);
		}
		return point;
	}

	@Override
	public String toString() {
		List<String> strings = new ArrayList<>();
		for (Isogeny<T> t : isogenies) {
			strings.add(t.toString());
		}
		return String.join(" o ", strings);
	}

	@Override
	public int compareTo(Isogeny<T> o) {
		if (o instanceof KernelIsogeny<?>) {
			KernelIsogeny<T> other = (KernelIsogeny<T>) o;
			if (isogenies.size() != other.isogenies.size()) {
				return isogenies.size() - other.isogenies.size();
			}
			Iterator<Isogeny<T>> it = isogenies.iterator();
			Iterator<Isogeny<T>> oIt = other.isogenies.iterator();
			while (it.hasNext()) {
				int cmp = it.next().compareTo(oIt.next());
				if (cmp != 0) {
					return cmp;
				}
			}
			return 0;
		}
		return getClass().getName().compareTo(o.getClass().getName());
	}

	@Override
	public EllipticCurve<T> getDomain() {
		return isogenies.getFirst().getDomain();
	}

	@Override
	public EllipticCurve<T> getRange() {
		return isogenies.getLast().getRange();
	}

	@Override
	public ProjectivePoint<T> evaluate(ProjectivePoint<T> point) {
		return evaluateLadder(point, isogenies);
	}

	@Override
	public BigInteger getDegree() {
		return degree.getValue();
	}

	@Override
	public List<ProjectivePoint<T>> kernelGenerators() {
		return kernelGenerators;
	}

	@Override
	public ProjectiveMorphism<T> asMorphism() {
		// TODO Auto-generated method stub
		return null;
	}

	private static class DualKernelIsogeny<T extends Element<T>> extends AbstractElement<Isogeny<T>>
			implements Isogeny<T> {
		private KernelIsogeny<T> dual;

		private DualKernelIsogeny(KernelIsogeny<T> dual) {
			this.dual = dual;
		}

		@Override
		public String toString() {
			return dual.toString() + "~";
		}

		@Override
		public int compareTo(Isogeny<T> o) {
			if (o instanceof DualKernelIsogeny<?>) {
				DualKernelIsogeny<T> other = (DualKernelIsogeny<T>) o;
				if (dual.isogenies.size() != other.dual.isogenies.size()) {
					return dual.isogenies.size() - other.dual.isogenies.size();
				}
				Iterator<Isogeny<T>> it = dual.isogenies.iterator();
				Iterator<Isogeny<T>> oIt = other.dual.isogenies.iterator();
				while (it.hasNext()) {
					int cmp = it.next().compareTo(oIt.next());
					if (cmp != 0) {
						return cmp;
					}
				}
				return 0;
			}
			return getClass().getName().compareTo(o.getClass().getName());
		}

		@Override
		public EllipticCurve<T> getDomain() {
			return dual.getRange();
		}

		@Override
		public EllipticCurve<T> getRange() {
			return dual.getDomain();
		}

		@Override
		public ProjectivePoint<T> evaluate(ProjectivePoint<T> point) {
			Iterator<Isogeny<T>> it = dual.isogenies.descendingIterator();
			while (it.hasNext()) {
				Isogeny<T> isogeny = it.next();
				point = isogeny.getDual().evaluate(point);
			}
			return point;
		}

		@Override
		public BigInteger getDegree() {
			return dual.getDegree();
		}

		@Override
		public Isogeny<T> getDual() {
			return dual;
		}

		@Override
		public ProjectiveMorphism<T> asMorphism() {
			// TODO Auto-generated method stub
			return null;
		}

		@Override
		public List<ProjectivePoint<T>> kernelGenerators() {
			// TODO Auto-generated method stub
			return null;
		}
	}

	@Override
	public Isogeny<T> getDual() {
		return new DualKernelIsogeny<>(this);
	}
}
