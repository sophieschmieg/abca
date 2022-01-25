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
import fields.interfaces.Element;
import varieties.projective.ProjectiveMorphism;
import varieties.projective.ProjectivePoint;

public class KernelIsogeny<T extends Element<T>> extends AbstractElement<Isogeny<T>> implements Isogeny<T> {
	private Deque<Isogeny<T>> isogenies;
	private BigInteger degree;

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
	}

	@Override
	public Isogeny<T> getDual() {
		return new DualKernelIsogeny<>(this);
	}
}
