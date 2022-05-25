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
import fields.interfaces.Module;
import fields.interfaces.Ring;
import fields.vectors.ExteriorProduct.WedgeProductSum;
import util.MiscAlgorithms;

public class ExteriorProduct<T extends Element<T>, S extends Element<S>>
		extends AbstractAlgebra<T, WedgeProductSum<T>> {
	public static class WedgeProduct implements Comparable<WedgeProduct> {
		private SortedSet<Integer> variables;
		private String[] variableNames;

		private WedgeProduct(Collection<Integer> variables, String[] variableNames) {
			this.variables = new TreeSet<>();
			this.variables.addAll(variables);
			this.variableNames = variableNames;
		}

		public SortedSet<Integer> getVariables() {
			return variables;
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
				build.append(variableNames[i - 1]);
			}
			return build.toString();
		}

		@Override
		public int compareTo(WedgeProduct o) {
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
		private Map<WedgeProduct, T> coefficients;
		private Ring<T> ring;

		private WedgeProductSum(Map<WedgeProduct, T> coefficients, Ring<T> ring) {
			this.coefficients = new TreeMap<>();
			for (WedgeProduct primitive : coefficients.keySet()) {
				T coefficient = coefficients.get(primitive);
				if (!coefficient.equals(ring.zero())) {
					this.coefficients.put(primitive, coefficient);
				}
			}
			this.ring = ring;
		}

		public Set<WedgeProduct> getWedgeProducts() {
			return coefficients.keySet();
		}

		@Override
		public String toString() {
			if (coefficients.isEmpty()) {
				return "0";
			}
			StringBuilder build = new StringBuilder();
			boolean first = true;
			for (WedgeProduct primitive : coefficients.keySet()) {
				if (first) {
					first = false;
				} else {
					build.append(" + ");
				}
				T coefficient = coefficient(primitive);
				boolean empty = primitive.variables.isEmpty();
				if (!coefficient.equals(ring.one())) {
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

		public T coefficient(WedgeProduct primitive) {
			return coefficients.getOrDefault(primitive, ring.zero());
		}

		@Override
		public int compareTo(WedgeProductSum<T> o) {
			Set<WedgeProduct> primitives = new TreeSet<>();
			primitives.addAll(coefficients.keySet());
			primitives.addAll(o.coefficients.keySet());
			for (WedgeProduct primitive : primitives) {
				T here = coefficient(primitive);
				T there = o.coefficient(primitive);
				int cmp = here.compareTo(there);
				if (cmp != 0) {
					return cmp;
				}
			}
			return 0;
		}
	}

	private Ring<T> ring;
	private Module<T, S> module;
	private String[] variableNames;
	private List<WedgeProductSum<T>> generators;
	private List<List<WedgeProductSum<T>>> gradedGenerators;
	private FreeModule<T> asFreeModule;
	private List<FreeModule<T>> asGradedFreeModule;

	private static <T extends Element<T>, S extends Element<S>> String[] defaultVariableNames(Module<T, S> module) {
		List<S> generators = module.getModuleGenerators();
		String[] result = new String[generators.size()];
		for (int i = 0; i < result.length; i++) {
			result[i] = generators.get(i).toString();
		}
		return result;
	}

	public ExteriorProduct(Module<T, S> module) {
		this(module, defaultVariableNames(module));
	}

	public ExteriorProduct(Module<T, S> module, String[] variableNames) {
		if (!module.isFree()) {
			throw new ArithmeticException("Not a free module!");
		}
		this.module = module;
		this.ring = module.getRing();
		this.variableNames = variableNames;
	}

	public FreeModule<T> asFreeModule() {
		if (asFreeModule == null) {
			this.asFreeModule = new FreeModule<>(module.getRing(),
					BigInteger.TWO.pow(module.getModuleGenerators().size()).intValueExact());
		}
		return asFreeModule;
	}
	
	public FreeModule<T> asGradedFreeModule(int degree) {
		if (degree < 0 || degree > variableNames.length) {
			throw new ArithmeticException("Undefined");
		}
		if (asGradedFreeModule == null) {
			asGradedFreeModule = new ArrayList<>();
		}
		while (asGradedFreeModule.size() <= degree) {
			asGradedFreeModule.add(null);
		}
		if (asGradedFreeModule.get(degree) == null) {
			asGradedFreeModule.set(degree, new FreeModule<>(ring, MiscAlgorithms.binomial(variableNames.length, degree)));
		}
		return asGradedFreeModule.get(degree);
	}

	@Override
	public Ring<T> getRing() {
		return ring;
	}

	@Override
	public WedgeProductSum<T> zero() {
		return new WedgeProductSum<>(Collections.emptyMap(), ring);
	}

	@Override
	public WedgeProductSum<T> add(WedgeProductSum<T> s1, WedgeProductSum<T> s2) {
		Set<WedgeProduct> primitives = new TreeSet<>();
		primitives.addAll(s1.coefficients.keySet());
		primitives.addAll(s2.coefficients.keySet());
		Map<WedgeProduct, T> result = new TreeMap<>();
		for (WedgeProduct primitive : primitives) {
			result.put(primitive, ring.add(s1.coefficient(primitive), s2.coefficient(primitive)));
		}
		return new WedgeProductSum<>(result, ring);
	}

	@Override
	public WedgeProductSum<T> negative(WedgeProductSum<T> s) {
		Map<WedgeProduct, T> result = new TreeMap<>();
		for (WedgeProduct primitive : s.coefficients.keySet()) {
			result.put(primitive, ring.negative(s.coefficient(primitive)));
		}
		return new WedgeProductSum<>(result, ring);
	}

	@Override
	public WedgeProductSum<T> scalarMultiply(T t, WedgeProductSum<T> s) {
		Map<WedgeProduct, T> result = new TreeMap<>();
		for (WedgeProduct primitive : s.coefficients.keySet()) {
			result.put(primitive, ring.multiply(t, s.coefficient(primitive)));
		}
		return new WedgeProductSum<>(result, ring);
	}

	@Override
	public boolean isFree() {
		return module.isFree();
	}

	@Override
	public Ideal<T> annihilator() {
		return module.getRing().getZeroIdeal();
	}
	
	@Override
	public List<Vector<T>> getSyzygies() {
		return Collections.emptyList();
	}

	private List<Vector<T>> asVectorList(List<WedgeProductSum<T>> s) {
		List<Vector<T>> result = new ArrayList<>();
		for (WedgeProductSum<T> t : s) {
			result.add(asVector(t));
		}
		return result;
	}

	@Override
	public boolean isLinearIndependent(List<WedgeProductSum<T>> s) {
		return asFreeModule().isLinearIndependent(asVectorList(s));
	}

	@Override
	public boolean isGeneratingModule(List<WedgeProductSum<T>> s) {
		return asFreeModule().isGeneratingModule(asVectorList(s));
	}

	@Override
	public List<Vector<T>> nonTrivialCombinations(List<WedgeProductSum<T>> s) {
		return asFreeModule().nonTrivialCombinations(asVectorList(s));
	}

	@Override
	public List<WedgeProductSum<T>> getModuleGenerators() {
		if (generators == null) {
			this.generators = new ArrayList<>();
			this.gradedGenerators = new ArrayList<>();
			int rank = module.getModuleGenerators().size();
			for (int i = 0; i <= rank; i++) {
				gradedGenerators.add(new ArrayList<>());
			}
			for (int i = 0; i < BigInteger.TWO.pow(rank).intValueExact(); i++) {
				Set<Integer> variables = new TreeSet<>();
				for (int j = 0; j < rank; j++) {
					if (((1 << j) & i) != 0) {
						variables.add(j + 1);
					}
				}
				WedgeProductSum<T> generator = new WedgeProductSum<>(
						Collections.singletonMap(new WedgeProduct(variables, variableNames), ring.one()), ring);
				generators.add(generator);
				gradedGenerators.get(variables.size()).add(generator);
			}
			for (int i = 0; i <= rank; i++) {
				Collections.sort(gradedGenerators.get(i));
			}
			Collections.sort(generators);
		}
		return generators;
	}

	public List<WedgeProductSum<T>> getGradedModuleGenerators(int degree) {
		getModuleGenerators();
		return gradedGenerators.get(degree);
	}

	@Override
	public Vector<T> asVector(WedgeProductSum<T> s) {
		List<T> result = new ArrayList<>();
		for (WedgeProductSum<T> generator : getModuleGenerators()) {
			result.add(s.coefficient(generator.coefficients.keySet().iterator().next()));
		}
		return new Vector<>(result);
	}

	public Vector<T> asGradedVector(WedgeProductSum<T> s, int degree) {
		List<T> result = new ArrayList<>();
		for (WedgeProductSum<T> generator : getGradedModuleGenerators(degree)) {
			result.add(s.coefficient(generator.coefficients.keySet().iterator().next()));
		}
		return new Vector<>(result);

	}

	public T asScalar(WedgeProductSum<T> s) {
		return s.coefficient(new WedgeProduct(Collections.emptySet(), variableNames));
	}

	@Override
	public Exactness exactness() {
		return module.exactness();
	}

	@Override
	public WedgeProductSum<T> getRandomElement() {
		return fromVector(asFreeModule().getRandomElement());
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
			Iterator<Vector<T>> it = asFreeModule().iterator();

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

	public WedgeProduct getWedgeProduct(Collection<Integer> variables) {
		return new WedgeProduct(variables, variableNames);
	}

	@Override
	public WedgeProductSum<T> getEmbedding(T t) {
		return new WedgeProductSum<>(
				Collections.singletonMap(new WedgeProduct(Collections.emptySet(), variableNames), t), ring);
	}

	public WedgeProductSum<T> getEmbedding(T t, WedgeProduct wedge) {
		return new WedgeProductSum<>(Collections.singletonMap(wedge, t), ring);
	}

	public WedgeProductSum<T> getModuleEmbedding(S t) {
		Vector<T> asVector = module.asVector(t);
		Map<WedgeProduct, T> result = new TreeMap<>();
		for (int i = 0; i < variableNames.length; i++) {
			result.put(new WedgeProduct(Collections.singleton(i + 1), variableNames), asVector.get(i + 1));
		}
		return new WedgeProductSum<>(result, ring);
	}

	public WedgeProductSum<T> getHigherModuleEmbedding(S t) {
		Vector<T> asVector = module.asVector(t);
		Map<WedgeProduct, T> result = new TreeMap<>();
		Set<Integer> allVariables = new TreeSet<>();
		for (int i = 0; i < variableNames.length; i++) {
			allVariables.add(i + 1);
		}
		for (int i = 0; i < variableNames.length; i++) {
			Set<Integer> variables = new TreeSet<>();
			variables.addAll(allVariables);
			variables.remove(i + 1);
			result.put(new WedgeProduct(variables, variableNames), asVector.get(i + 1));
		}
		return new WedgeProductSum<>(result, ring);
	}

	private WedgeProductSum<T> wedge(WedgeProduct t1, WedgeProduct t2) {
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
				Collections.singletonMap(new WedgeProduct(common, variableNames), ring.getInteger(signum)), ring);
	}

	private WedgeProductSum<T> koszulDerivative(WedgeProduct t, List<T> generators) {
		int sign = 1;
		Map<WedgeProduct, T> result = new TreeMap<>();
		for (int variable : t.variables) {
			Set<Integer> removed = new TreeSet<>();
			removed.addAll(t.variables);
			removed.remove(variable);
			result.put(new WedgeProduct(removed, variableNames), ring.multiply(sign, generators.get(variable - 1)));
			sign *= -1;
		}
		return new WedgeProductSum<>(result, ring);
	}

	public WedgeProductSum<T> koszulDerivative(WedgeProductSum<T> t, List<T> generators) {
		WedgeProductSum<T> result = zero();
		for (WedgeProduct primitive : t.coefficients.keySet()) {
			result = add(scalarMultiply(t.coefficient(primitive), koszulDerivative(primitive, generators)), result);
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
		for (S generator : module.getModuleGenerators()) {
			result.add(getModuleEmbedding(generator));
		}
		return result;
	}

	@Override
	public WedgeProductSum<T> one() {
		return getEmbedding(ring.one());
	}

	@Override
	public BigInteger characteristic() {
		return ring.characteristic();
	}

	@Override
	public WedgeProductSum<T> multiply(WedgeProductSum<T> t1, WedgeProductSum<T> t2) {
		WedgeProductSum<T> result = zero();
		for (WedgeProduct primitive1 : t1.coefficients.keySet()) {
			for (WedgeProduct primitive2 : t2.coefficients.keySet()) {
				result = add(scalarMultiply(ring.multiply(t1.coefficient(primitive1), t2.coefficient(primitive2)),
						wedge(primitive1, primitive2)), result);
			}
		}
		return result;
	}

	@Override
	public boolean isUnit(WedgeProductSum<T> t) {
		if (t.coefficients.size() != 1) {
			return false;
		}
		return ring.isUnit(t.coefficient(new WedgeProduct(Collections.emptySet(), variableNames)));
	}

	@Override
	public WedgeProductSum<T> inverse(WedgeProductSum<T> t) {
		if (!isUnit(t)) {
			throw new ArithmeticException("Not a unit!");
		}
		return getEmbedding(ring.inverse(t.coefficient(new WedgeProduct(Collections.emptySet(), variableNames))));
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
		return ring.isIrreducible();
	}

	@Override
	public boolean isZeroDivisor(WedgeProductSum<T> t) {
		for (WedgeProduct primitive : t.coefficients.keySet()) {
			if (!primitive.variables.isEmpty()) {
				return true;
			}
			if (ring.isZeroDivisor(t.coefficient(primitive))) {
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
		WedgeProduct leading = t.coefficients.keySet().iterator().next();
		T leadingCoefficient = t.coefficient(leading);
		return scalarMultiply(ring.inverse(ring.projectToUnit(leadingCoefficient)), t);
	}

	@Override
	public Iterable<WedgeProductSum<T>> getUnits() {
		return new Iterable<>() {

			@Override
			public Iterator<WedgeProductSum<T>> iterator() {
				return new Iterator<>() {
					private Iterator<T> it = ring.getUnits().iterator();

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
