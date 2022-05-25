package fields.vectors;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import fields.exceptions.InfinityException;
import fields.helper.AbstractAlgebra;
import fields.helper.AbstractElement;
import fields.interfaces.Element;
import fields.interfaces.Ideal;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.polynomials.PolynomialExteriorProduct.WedgeProductSum;
import fields.vectors.FreeModule;
import fields.vectors.Vector;

public class PolynomialExteriorProduct<T extends Element<T>>
		extends AbstractAlgebra<Polynomial<T>, WedgeProductSum<T>> {
	public static class WedgeProduct<T extends Element<T>> implements Comparable<WedgeProduct<T>> {
		private SortedSet<Integer> variables;
		private String[] variableNames;

		private WedgeProduct(Collection<Integer> variables, String[] variableNames) {
			this.variables = new TreeSet<>();
			this.variables.addAll(variables);
			this.variableNames = variableNames;
		}

		@Override
		public String toString() {
			StringBuilder build = new StringBuilder();
			boolean first = true;
			for (int i : variables) {
				if (first) {
					first = false;
				} else {
					build.append(" ^ ");
				}
				build.append("d");
				build.append(variableNames[i - 1]);
			}
			return build.toString();
		}

		@Override
		public int compareTo(WedgeProduct<T> o) {
			if (variables.size() != o.variables.size()) {
				return o.variables.size() - variables.size();
			}
			for (int i = 0; i < variableNames.length; i++) {
				boolean hasI = variables.contains(i + 1);
				boolean otherHasI = o.variables.contains(i + 1);
				if (hasI && !otherHasI) {
					return 1;
				}
				if (!hasI && otherHasI) {
					return -1;
				}
			}
			return 0;
		}
	}

	public static class WedgeProductSum<T extends Element<T>> extends AbstractElement<WedgeProductSum<T>> {
		private Map<WedgeProduct<T>, Polynomial<T>> coefficients;
		private PolynomialRing<T> polynomialRing;

		private WedgeProductSum(Map<WedgeProduct<T>, Polynomial<T>> coefficients, PolynomialRing<T> polynomialRing) {
			this.coefficients = new TreeMap<>();
			for (WedgeProduct<T> primitive : coefficients.keySet()) {
				Polynomial<T> coefficient = coefficients.get(primitive);
				if (!coefficient.equals(polynomialRing.zero())) {
					this.coefficients.put(primitive, coefficient);
				}
			}
			this.polynomialRing = polynomialRing;
		}

		@Override
		public String toString() {
			if (coefficients.isEmpty()) {
				return "0";
			}
			StringBuilder build = new StringBuilder();
			boolean first = true;
			for (WedgeProduct<T> primitive : coefficients.keySet()) {
				if (first) {
					first = false;
				} else {
					build.append(" + ");
				}
				Polynomial<T> coefficient = coefficient(primitive);
				boolean empty = primitive.variables.isEmpty();
				if (!coefficient.equals(polynomialRing.one())) {
					String asString = coefficient.toString();
					boolean space = asString.contains(" ") && !empty;
					build.append((space ? "(" : "") + asString + (space ? ")" : "") + (empty ? "" : "*"));
				} else if (empty) {
					build.append("1");
				}
				build.append(primitive);
			}
			return build.toString();
		}

		public Polynomial<T> coefficient(WedgeProduct<T> primitive) {
			return coefficients.getOrDefault(primitive, polynomialRing.zero());
		}

		@Override
		public int compareTo(WedgeProductSum<T> o) {
			Set<WedgeProduct<T>> primitives = new TreeSet<>();
			primitives.addAll(coefficients.keySet());
			primitives.addAll(o.coefficients.keySet());
			for (WedgeProduct<T> primitive : primitives) {
				Polynomial<T> here = coefficient(primitive);
				Polynomial<T> there = o.coefficient(primitive);
				int cmp = here.compareTo(there);
				if (cmp != 0) {
					return cmp;
				}
			}
			return 0;
		}
	}

	private PolynomialRing<T> polynomialRing;
	private List<WedgeProductSum<T>> generators;
	private List<List<WedgeProductSum<T>>> gradedGenerators;
	private FreeModule<Polynomial<T>> asFreeModule;

	public PolynomialExteriorProduct(PolynomialRing<T> polynomialRing) {
		this.polynomialRing = polynomialRing;
		this.asFreeModule = new FreeModule<>(polynomialRing,
				BigInteger.TWO.pow(polynomialRing.numberOfVariables()).intValueExact());
		this.generators = new ArrayList<>();
		this.gradedGenerators = new ArrayList<>();
		for (int i = 0; i <= polynomialRing.numberOfVariables(); i++) {
			gradedGenerators.add(new ArrayList<>());
		}
		for (int i = 0; i < BigInteger.TWO.pow(polynomialRing.numberOfVariables()).intValueExact(); i++) {
			Set<Integer> variables = new TreeSet<>();
			for (int j = 0; j < polynomialRing.numberOfVariables(); j++) {
				if (((1 << j) & i) != 0) {
					variables.add(j + 1);
				}
			}
			WedgeProductSum<T> generator = new WedgeProductSum<>(
					Collections.singletonMap(new WedgeProduct<>(variables, polynomialRing.getVariableNames()),
							polynomialRing.one()),
					polynomialRing);
			generators.add(generator);
			gradedGenerators.get(variables.size()).add(generator);
		}
		for (int i = 0; i <= polynomialRing.numberOfVariables(); i++) {
			Collections.sort(gradedGenerators.get(i));
		}
		Collections.sort(generators);
	}

	@Override
	public PolynomialRing<T> getRing() {
		return polynomialRing;
	}

	@Override
	public WedgeProductSum<T> zero() {
		return new WedgeProductSum<>(Collections.emptyMap(), polynomialRing);
	}

	@Override
	public WedgeProductSum<T> add(WedgeProductSum<T> s1, WedgeProductSum<T> s2) {
		Set<WedgeProduct<T>> primitives = new TreeSet<>();
		primitives.addAll(s1.coefficients.keySet());
		primitives.addAll(s2.coefficients.keySet());
		Map<WedgeProduct<T>, Polynomial<T>> result = new TreeMap<>();
		for (WedgeProduct<T> primitive : primitives) {
			result.put(primitive, polynomialRing.add(s1.coefficient(primitive), s2.coefficient(primitive)));
		}
		return new WedgeProductSum<>(result, polynomialRing);
	}

	@Override
	public WedgeProductSum<T> negative(WedgeProductSum<T> s) {
		Map<WedgeProduct<T>, Polynomial<T>> result = new TreeMap<>();
		for (WedgeProduct<T> primitive : s.coefficients.keySet()) {
			result.put(primitive, polynomialRing.negative(s.coefficient(primitive)));
		}
		return new WedgeProductSum<>(result, polynomialRing);
	}

	@Override
	public WedgeProductSum<T> scalarMultiply(Polynomial<T> t, WedgeProductSum<T> s) {
		Map<WedgeProduct<T>, Polynomial<T>> result = new TreeMap<>();
		for (WedgeProduct<T> primitive : s.coefficients.keySet()) {
			result.put(primitive, polynomialRing.multiply(t, s.coefficient(primitive)));
		}
		return new WedgeProductSum<>(result, polynomialRing);
	}

	@Override
	public boolean isFree() {
		return true;
	}

	@Override
	public PolynomialIdeal<T> annihilator() {
		return polynomialRing.getZeroIdeal();
	}

	private List<Vector<Polynomial<T>>> asVectorList(List<WedgeProductSum<T>> s) {
		List<Vector<Polynomial<T>>> result = new ArrayList<>();
		for (WedgeProductSum<T> t : s) {
			result.add(asVector(t));
		}
		return result;
	}

	@Override
	public boolean isLinearIndependent(List<WedgeProductSum<T>> s) {
		return asFreeModule.isLinearIndependent(asVectorList(s));
	}

	@Override
	public boolean isGeneratingModule(List<WedgeProductSum<T>> s) {
		return asFreeModule.isGeneratingModule(asVectorList(s));
	}

	@Override
	public List<List<Polynomial<T>>> nonTrivialCombinations(List<WedgeProductSum<T>> s) {
		return asFreeModule.nonTrivialCombinations(asVectorList(s));
	}

	@Override
	public List<WedgeProductSum<T>> getModuleGenerators() {
		return generators;
	}

	public List<WedgeProductSum<T>> getGradedModuleGenerators(int degree) {
		return gradedGenerators.get(degree);
	}

	@Override
	public Vector<Polynomial<T>> asVector(WedgeProductSum<T> s) {
		List<Polynomial<T>> result = new ArrayList<>();
		for (WedgeProductSum<T> generator : generators) {
			result.add(s.coefficient(generator.coefficients.keySet().iterator().next()));
		}
		return new Vector<>(result);
	}

	@Override
	public Exactness exactness() {
		return polynomialRing.exactness();
	}

	@Override
	public WedgeProductSum<T> getRandomElement() {
		return fromVector(asFreeModule.getRandomElement());
	}

	@Override
	public boolean isFinite() {
		return false;
	}

	@Override
	public BigInteger getNumberOfElements() throws InfinityException {
		throw new InfinityException();
	}

	@Override
	public Iterator<WedgeProductSum<T>> iterator() {
		return new Iterator<>() {
			Iterator<Vector<Polynomial<T>>> it = asFreeModule.iterator();

			@Override
			public boolean hasNext() {
				return it.hasNext();
			}

			@Override
			public WedgeProductSum<T> next() {
				return fromVector(it.next());
			}
		};
	}

	@Override
	public WedgeProductSum<T> getEmbedding(Polynomial<T> t) {
		return new WedgeProductSum<>(Collections.singletonMap(
				new WedgeProduct<>(Collections.emptySet(), polynomialRing.getVariableNames()), t), polynomialRing);
	}

	private WedgeProductSum<T> wedge(WedgeProduct<T> t1, WedgeProduct<T> t2) {
		Set<Integer> intersection = new TreeSet<>();
		intersection.addAll(t1.variables);
		intersection.retainAll(t2.variables);
		if (!intersection.isEmpty()) {
			return zero();
		}
		int signum = 1;
		List<Integer> common = new ArrayList<>();
		common.addAll(t1.variables);
		common.addAll(t2.variables);
		for (int i = 0; i < common.size(); i++) {
			for (int j = 0; j < i; j++) {
				if (common.get(j) > common.get(i)) {
					signum *= -1;
				}
			}
		}
		return new WedgeProductSum<>(
				Collections.singletonMap(new WedgeProduct<>(common, polynomialRing.getVariableNames()),
						polynomialRing.getInteger(signum)),
				polynomialRing);
	}

	private WedgeProductSum<T> koszulDerivative(WedgeProduct<T> t, List<Polynomial<T>> generators) {
		int sign = 1;
		Map<WedgeProduct<T>, Polynomial<T>> result = new TreeMap<>();
		for (int variable : t.variables) {
			Set<Integer> removed = new TreeSet<>();
			removed.addAll(t.variables);
			removed.remove(variable);
			result.put(new WedgeProduct<>(removed, polynomialRing.getVariableNames()),
					polynomialRing.multiply(sign, generators.get(variable - 1)));
			sign *= -1;
		}
		return new WedgeProductSum<>(result, polynomialRing);
	}

	public WedgeProductSum<T> koszulDerivative(WedgeProductSum<T> t, List<Polynomial<T>> generators) {
		WedgeProductSum<T> result = zero();
		for (WedgeProduct<T> primitive : t.coefficients.keySet()) {
			result = add(scalarMultiply(t.coefficient(primitive), koszulDerivative(primitive, generators)), result);
		}
		return result;
	}

	public WedgeProductSum<T> derivative(Polynomial<T> t) {
		return derivative(getEmbedding(t));
	}

	public WedgeProductSum<T> derivative(WedgeProductSum<T> t) {
		WedgeProductSum<T> result = zero();
		for (WedgeProduct<T> primitive : t.coefficients.keySet()) {
			Polynomial<T> coefficient = t.coefficient(primitive);
			for (int i = 0; i < polynomialRing.numberOfVariables(); i++) {
				if (primitive.variables.contains(i + 1)) {
					continue;
				}
				result = add(scalarMultiply(polynomialRing.derivative(coefficient, i + 1),
						wedge(new WedgeProduct<>(Collections.singleton(i + 1), polynomialRing.getVariableNames()),
								primitive)),
						result);
			}
		}
		return result;
	}

	@Override
	public boolean isGeneratingAlgebra(List<WedgeProductSum<T>> s) {
		throw new ArithmeticException("Not implemented");
	}

	@Override
	public List<WedgeProductSum<T>> getAlgebraGenerators() {
		List<WedgeProductSum<T>> result = new ArrayList<>();
		for (int i = 0; i < polynomialRing.numberOfVariables(); i++) {
			result.add(derivative(polynomialRing.getVar(i + 1)));
		}
		return result;
	}

	@Override
	public WedgeProductSum<T> one() {
		return getEmbedding(polynomialRing.one());
	}

	@Override
	public BigInteger characteristic() {
		return polynomialRing.characteristic();
	}

	@Override
	public WedgeProductSum<T> multiply(WedgeProductSum<T> t1, WedgeProductSum<T> t2) {
		WedgeProductSum<T> result = zero();
		for (WedgeProduct<T> primitive1 : t1.coefficients.keySet()) {
			for (WedgeProduct<T> primitive2 : t2.coefficients.keySet()) {
				result = add(
						scalarMultiply(polynomialRing.multiply(t1.coefficient(primitive1), t2.coefficient(primitive2)),
								wedge(primitive1, primitive2)),
						result);
			}
		}
		return result;
	}

	@Override
	public boolean isUnit(WedgeProductSum<T> t) {
		if (t.coefficients.size() != 1) {
			return false;
		}
		return polynomialRing
				.isUnit(t.coefficient(new WedgeProduct<>(Collections.emptySet(), polynomialRing.getVariableNames())));
	}

	@Override
	public WedgeProductSum<T> inverse(WedgeProductSum<T> t) {
		if (!isUnit(t)) {
			throw new ArithmeticException("Not a unit!");
		}
		return getEmbedding(polynomialRing
				.inverse(t.coefficient(new WedgeProduct<>(Collections.emptySet(), polynomialRing.getVariableNames()))));
	}

	@Override
	public boolean isCommutative() {
		return false;
	}

	@Override
	public boolean isIntegral() {
		return false;
	}

	@Override
	public boolean isReduced() {
		return false;
	}

	@Override
	public boolean isIrreducible() {
		return polynomialRing.isIrreducible();
	}

	@Override
	public boolean isZeroDivisor(WedgeProductSum<T> t) {
		for (WedgeProduct<T> primitive : t.coefficients.keySet()) {
			if (!primitive.variables.isEmpty()) {
				return true;
			}
			if (polynomialRing.isZeroDivisor(t.coefficient(primitive))) {
				return true;
			}
		}
		return false;
	}

	@Override
	public boolean isEuclidean() {
		return false;
	}

	@Override
	public boolean isUniqueFactorizationDomain() {
		return false;
	}

	@Override
	public FactorizationResult<WedgeProductSum<T>, WedgeProductSum<T>> uniqueFactorization(WedgeProductSum<T> t) {
		throw new ArithmeticException("Not a UFD!");
	}

	@Override
	public boolean isPrincipalIdealDomain() {
		return false;
	}

	@Override
	public boolean isDedekindDomain() {
		return false;
	}

	@Override
	public boolean isDivisible(WedgeProductSum<T> dividend, WedgeProductSum<T> divisor) {
		throw new ArithmeticException("No.");
	}

	@Override
	public QuotientAndRemainderResult<WedgeProductSum<T>> quotientAndRemainder(WedgeProductSum<T> dividend,
			WedgeProductSum<T> divisor) {
		throw new ArithmeticException("No.");
	}

	@Override
	public BigInteger euclidMeasure(WedgeProductSum<T> t) {
		throw new ArithmeticException("No.");
	}

	@Override
	public WedgeProductSum<T> projectToUnit(WedgeProductSum<T> t) {
		if (t.equals(zero())) {
			return one();
		}
		WedgeProduct<T> leading = t.coefficients.keySet().iterator().next();
		Polynomial<T> leadingPolynomial = t.coefficient(leading);
		return scalarMultiply(polynomialRing.inverse(polynomialRing.projectToUnit(leadingPolynomial)), t);
	}

	@Override
	public Iterable<WedgeProductSum<T>> getUnits() {
		return new Iterable<>() {

			@Override
			public Iterator<WedgeProductSum<T>> iterator() {
				return new Iterator<>() {
					private Iterator<Polynomial<T>> it = polynomialRing.getUnits().iterator();

					@Override
					public boolean hasNext() {
						return it.hasNext();
					}

					@Override
					public WedgeProductSum<T> next() {
						return getEmbedding(it.next());
					}
				};
			}
		};
	}

	@Override
	public int krullDimension() {
		throw new ArithmeticException("No.");
	}

	@Override
	public IdealResult<WedgeProductSum<T>, ?> getIdealWithTransforms(List<WedgeProductSum<T>> generators) {
		throw new ArithmeticException("No.");
	}

	@Override
	public Ideal<WedgeProductSum<T>> intersect(Ideal<WedgeProductSum<T>> t1, Ideal<WedgeProductSum<T>> t2) {
		throw new ArithmeticException("No.");
	}

	@Override
	public Ideal<WedgeProductSum<T>> radical(Ideal<WedgeProductSum<T>> t) {
		throw new ArithmeticException("No.");
	}
}
