package varieties.curves;

import java.math.BigInteger;
import java.util.Collections;
import java.util.Deque;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import fields.Element;
import varieties.ProjectivePoint;

public class KernelIsogeny<T extends Element> implements Isogeny<T> {
	private Deque<Isogeny<T>> isogenies;
	private BigInteger degree;

	public KernelIsogeny(EllipticCurve<T> domain, ProjectivePoint<T> kernel) {
		this(domain, Collections.singletonMap(kernel, domain.getNumberOfElements()));
	}

	public KernelIsogeny(EllipticCurve<T> domain, ProjectivePoint<T> kernel, BigInteger order) {
		this(domain, Collections.singletonMap(kernel, order));
	}

	private static <T extends Element> Map<ProjectivePoint<T>, BigInteger> fillMap(BigInteger order,
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
		this.isogenies = new LinkedList<>();
		EllipticCurve<T> currentRange = domain;
		for (ProjectivePoint<T> kernelGenerator : kernel.keySet()) {
			BigInteger order = kernel.get(kernelGenerator);
			if (!domain.hasRationalPoint(kernelGenerator)) {
				throw new ArithmeticException("Point " + kernelGenerator + " not on curve!");
			}
			if (!domain.multiply(order, kernelGenerator).equals(domain.neutral())) {
				throw new ArithmeticException("Order hint " + order + " for " + kernelGenerator + " was wrong!");
			}
			BigInteger prime = BigInteger.TWO;
			Map<BigInteger, Integer> orderPrimeDecomposition = new TreeMap<>();
			BigInteger reducedOrder = order;
			while (!reducedOrder.equals(BigInteger.ONE)) {
				if (!reducedOrder.mod(prime).equals(BigInteger.ZERO)) {
					prime = prime.nextProbablePrime();
					continue;
				}
				reducedOrder = order.divide(prime);
				if (domain.multiply(reducedOrder, kernelGenerator).equals(domain.neutral())) {
					order = reducedOrder;
					continue;
				}
				int i = 1;
				while (reducedOrder.mod(prime).equals(BigInteger.ZERO)) {
					i++;
					reducedOrder = reducedOrder.divide(prime);
				}
				orderPrimeDecomposition.put(prime, i);
				prime = prime.nextProbablePrime();
			}
			ProjectivePoint<T> kernelGeneratorImage = this.evaluate(kernelGenerator);
			for (BigInteger orderPrime : orderPrimeDecomposition.keySet()) {
				int exponent = orderPrimeDecomposition.get(orderPrime);
				BigInteger primePower = orderPrime.pow(exponent);
				while (!primePower.equals(BigInteger.ONE)) {
					primePower = primePower.divide(orderPrime);
					ProjectivePoint<T> powerKernelPoint = currentRange.multiply(primePower, kernelGeneratorImage);
					Isogeny<T> nextIsogeny = new KernelPointIsogeny<T>(currentRange, powerKernelPoint,
							orderPrime.intValueExact());
					this.isogenies.add(nextIsogeny);
					currentRange = nextIsogeny.getRange();
					kernelGeneratorImage = nextIsogeny.evaluate(kernelGeneratorImage);
				}
			}
		}
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
		for (Isogeny<T> isogeny : isogenies) {
			point = isogeny.evaluate(point);
		}
		return point;
	}

	@Override
	public BigInteger getDegree() {
		return degree;
	}

	@Override
	public Isogeny<T> getDual() {
		return new Isogeny<T>() {

			@Override
			public EllipticCurve<T> getDomain() {
				return KernelIsogeny.this.getRange();
			}

			@Override
			public EllipticCurve<T> getRange() {
				return KernelIsogeny.this.getDomain();
			}

			@Override
			public ProjectivePoint<T> evaluate(ProjectivePoint<T> point) {
				Iterator<Isogeny<T>> it = isogenies.descendingIterator();
				while (it.hasNext()) {
					Isogeny<T> isogeny = it.next();
					point = isogeny.evaluate(point);
				}
				return point;
			}

			@Override
			public BigInteger getDegree() {
				return KernelIsogeny.this.getDegree();
			}

			@Override
			public Isogeny<T> getDual() {
				return KernelIsogeny.this;
			}
		};
	}

}
