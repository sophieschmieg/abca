package fields.polynomials;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import fields.exceptions.InfinityException;
import fields.helper.AbstractAlgebra;
import fields.interfaces.Element;
import fields.interfaces.Field;
import fields.interfaces.Ideal;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.Ring;
import varieties.AffinePoint;

public abstract class AbstractPolynomialRing<T extends Element<T>> extends AbstractAlgebra<T, Polynomial<T>>
		implements PolynomialRing<T> {

	AbstractPolynomialRing(boolean generateUnivariatePolynomialRing) {
		super(generateUnivariatePolynomialRing);
	}

	public static <T extends Element<T>> AbstractPolynomialRing<T> getPolynomialRing(Ring<T> ring, int numvars,
			Comparator<Monomial> comparator) {
		if (numvars == 1) {
			return new UnivariatePolynomialRing<T>(ring);
		}
		return new MultivariatePolynomialRing<T>(ring, numvars, comparator);
	}

	@Override
	public Polynomial<T> projectToUnit(Polynomial<T> t) {
		return getEmbedding(getRing().projectToUnit(t.leadingCoefficient()));
	}

	@Override
	public Polynomial<T> upToUnit(Polynomial<T> t) {
		return divide(t, getRing().projectToUnit(t.leadingCoefficient()));
	}

	@Override
	public Polynomial<T> multiply(T t, Monomial m, Polynomial<T> p) {
		return multiply(getEmbedding(t, m.exponents()), p);
	}

	@Override
	public MultivariatePolynomial<T> homogenize(Polynomial<T> t) {
		Map<Monomial, T> coeff = new TreeMap<>();
		for (Monomial m : t.monomials()) {
			coeff.put(m.homogenizeMonomial(t.degree()), t.coefficient(m));
		}
		return new MultivariatePolynomial<>(coeff,
				new MultivariatePolynomialRing<>(getRing(), numberOfVariables() + 1, getComparator()));
	}

	@Override
	public Polynomial<T> dehomogenize(Polynomial<T> t, int coord) {
		if (coord < 1 || coord > this.numberOfVariables() || !isHomogeneous(t))
			throw new ArithmeticException("Not possible");
		List<T> values = new ArrayList<T>();
		for (int i = 0; i < this.numberOfVariables(); i++)
			if (i == coord - 1)
				values.add(this.getRing().one());
			else
				values.add(null);
		return this.partiallyEvaluate(t, values);
	}

	@SafeVarargs
	public final Polynomial<T> getLinear(T... coeffs) {
		return this.getLinear(Arrays.asList(coeffs));
	}

	@SafeVarargs
	public final T evaluate(Polynomial<T> t, T... values) {
		return this.evaluate(t, Arrays.asList(values));
	}

	@Override
	public int krullDimension() {
		return getRing().krullDimension() + 1;
	}

	@Override
	public List<Polynomial<T>> quotientAndRemainder(Polynomial<T> dividend, Polynomial<T> divisor) {
		return generalQuotientAndRemainder(dividend, Collections.singletonList(divisor));
	}

	@Override
	public Polynomial<T> reduce(Polynomial<T> polynomial, List<Polynomial<T>> basis) {
		if (basis.isEmpty()) {
			return polynomial;
		}
		return generalQuotientAndRemainder(polynomial, basis).get(basis.size());
	}

	@Override
	public Ideal<Polynomial<T>> getIdeal(List<Polynomial<T>> generators) {
		return new PolynomialIdeal<T>(this, generators);
	}

	@Override
	public PolynomialIdeal<T> intersect(Ideal<Polynomial<T>> t1, Ideal<Polynomial<T>> t2) {
		List<Polynomial<T>> generators = new ArrayList<>();
		PolynomialRing<T> ring = this.addVariableWithElimination(1);
		Polynomial<T> t = ring.getVar(1);
		Polynomial<T> omt = ring.subtract(ring.one(), t);
		for (Polynomial<T> f : t1.generators()) {
			generators.add(ring.multiply(t, ring.getEmbeddingWithElimination(f, 1)));
		}
		for (Polynomial<T> g : t2.generators()) {
			generators.add(ring.multiply(omt, ring.getEmbeddingWithElimination(g, 1)));
		}
		Ideal<Polynomial<T>> intersectionIdeal = ring.getIdeal(generators);
		List<Polynomial<T>> intersectionGenerators = new ArrayList<Polynomial<T>>();
		for (Polynomial<T> b : intersectionIdeal.generators()) {
			if (b.leadingMonomial().exponents()[0] == 0) {
				intersectionGenerators.add(getEmbeddingWithElimination(b, -1));
			}
		}
		return new PolynomialIdeal<>(this, intersectionGenerators);

	}

	@Override
	public List<Polynomial<T>> buchberger(List<Polynomial<T>> generators) {
		List<Polynomial<T>> groebner = new ArrayList<Polynomial<T>>();
		for (Polynomial<T> basispolynomial : generators) {
			Polynomial<T> reduced = this.reduce(basispolynomial, groebner);
			if (!reduced.equals(this.zero())) {
				groebner.add(upToUnit(reduced));
			}
		}
		while (true) {
			List<Polynomial<T>> newGroebner = new ArrayList<>();
			for (int i = 1; i < groebner.size(); i++) {
				for (int j = 0; j < i; j++) {
					Polynomial<T> pi = groebner.get(i);
					Polynomial<T> pj = groebner.get(j);
					Monomial li = pi.leadingMonomial();
					Monomial lj = pj.leadingMonomial();
					T ci = pi.leadingCoefficient();
					T cj = pj.leadingCoefficient();
					T coefficient = getRing().lcm(pi.leadingCoefficient(), pj.leadingCoefficient());
					Monomial lcm = li.lcm(lj);
					li = lcm.divide(li);
					lj = lcm.divide(lj);
					ci = getRing().divide(coefficient, ci);
					cj = getRing().divide(coefficient, cj);
					Polynomial<T> pli = getEmbedding(ci, li.exponents());
					Polynomial<T> plj = getEmbedding(cj, lj.exponents());
					Polynomial<T> newGenerator = subtract(multiply(pli, pi), this.multiply(plj, pj));
					newGenerator = reduce(newGenerator, groebner);
					newGenerator = reduce(newGenerator, newGroebner);
					if (!newGenerator.equals(this.zero())) {
						newGroebner.add(upToUnit(newGenerator));
					}
					if (!getRing().isEuclidean() || getRing().krullDimension() == 0) {
						continue;
					}
					Monomial libyj = li.divide(lj);
					if (libyj != null) {
						List<T> extendedEuclid = getRing().extendedEuclidean(ci, cj);
						if (!extendedEuclid.get(0).equals(ci)) {
							newGroebner.add(subtract(multiply(extendedEuclid.get(1), pi),
									multiply(extendedEuclid.get(2), libyj, pj)));
						}
					}
					Monomial ljbyi = lj.divide(li);
					if (ljbyi != null) {
						List<T> extendedEuclid = getRing().extendedEuclidean(cj, ci);
						if (!extendedEuclid.get(0).equals(cj)) {
							newGroebner.add(subtract(multiply(extendedEuclid.get(1), pj),
									multiply(extendedEuclid.get(2), ljbyi, pi)));
						}
					}
				}
			}
			newGroebner = reduceBasis(newGroebner);
			newGroebner.addAll(groebner);
			if (groebner.equals(newGroebner)) {
				break;
			}
			groebner = newGroebner;
		}
		return reduceBasis(groebner);
	}

	@Override
	public List<Polynomial<T>> reduceBasis(List<Polynomial<T>> generators) {
		List<Polynomial<T>> list = new ArrayList<Polynomial<T>>();
		list.addAll(generators);
		SortedSet<Polynomial<T>> reduced = new TreeSet<Polynomial<T>>(Collections.reverseOrder());
		for (int i = 0; i < generators.size(); i++) {
			list.remove(i);
			Polynomial<T> reduce = this.reduce(generators.get(i), list);
			list.add(i, generators.get(i));
			if (!reduce.equals(this.zero())) {
				reduced.add(reduce);
			}
		}
		list.clear();
		list.addAll(reduced);
		return list;
	}

	@Override
	public boolean isGeneratingAlgebra(List<Polynomial<T>> s) {
		Ideal<Polynomial<T>> ideal = getIdeal(s);
		List<Polynomial<T>> variables = new ArrayList<>();
		for (int i = 0; i < numberOfVariables(); i++) {
			variables.add(getVar(i + 1));
		}
		Ideal<Polynomial<T>> all = getIdeal(variables);
		return ideal.equalsIdeal(all);
	}

	@Override
	public String toString() {
		if (numberOfVariables() == 1) {
			return getRing().toString() + "[X]";
		} else if (numberOfVariables() == 2) {
			return getRing().toString() + "[X,Y]";
		} else if (numberOfVariables() == 3) {
			return getRing().toString() + "[X,Y,Z]";
		}
		StringBuffer buf = new StringBuffer();
		buf.append(getRing().toString() + "[");
		boolean first = true;
		for (int i = 0; i < numberOfVariables(); i++) {
			if (first) {
				first = false;
			} else {
				buf.append(",");
			}
			buf.append("X_" + i);
		}
		buf.append("]");
		return buf.toString();
	}

	@Override
	public List<AffinePoint<T>> solve(List<Polynomial<T>> polynomials) {
		if (!(getRing() instanceof Field<?>)) {
			throw new UnsupportedOperationException();
		}
		Field<T> field = (Field<T>) getRing();
		if (polynomials.isEmpty()) {
			throw new InfinityException();
		}
		if (numberOfVariables() == 1) {
			Set<T> roots = new TreeSet<>();
			roots.addAll(field.roots(polynomials.get(0)));
			for (int i = 1; i < polynomials.size(); i++) {
				roots.retainAll(field.roots(polynomials.get(i)));
			}
			List<AffinePoint<T>> result = new ArrayList<>();
			for (T root : roots) {
				result.add(new AffinePoint<T>(field, Collections.singletonList(root)));
			}
			return result;
		}
		if (getComparator() != Monomial.REVLEX) {
			PolynomialRing<T> revlexRing = getPolynomialRing(getRing(), numberOfVariables(), Monomial.REVLEX);
			List<Polynomial<T>> embeded = new ArrayList<>();
			int[] map = new int[numberOfVariables()];
			for (int i = 0; i < numberOfVariables(); i++) {
				map[i] = i;
			}
			for (Polynomial<T> polynomial : polynomials) {
				embeded.add(revlexRing.getEmbedding(polynomial, map));
			}
			return revlexRing.solve(embeded);
		}
		Ideal<Polynomial<T>> eliminationIdeal = getIdeal(polynomials);
		List<Polynomial<T>> generators = eliminationIdeal.generators();
		int variable = 0;
		int highest = 0;
		List<AffinePoint<T>> prevPartialSolutions = new ArrayList<>();
		List<AffinePoint<T>> partialSolutions = new ArrayList<>();
		partialSolutions.add(new AffinePoint<>(field, Collections.emptyList()));
		int[] map = new int[numberOfVariables()];
		for (int i = 0; i < map.length; i++) {
			map[i] = -1;
		}
		for (int i = generators.size() - 1; i >= 0; i--) {
			Polynomial<T> p = generators.get(i);
			highest = highestVariable(p);
			if (highest <= variable) {
				prevPartialSolutions.clear();
				prevPartialSolutions.addAll(partialSolutions);
				partialSolutions.clear();
				for (AffinePoint<T> prevSolution : prevPartialSolutions) {
					List<T> values = new ArrayList<>();
					values.addAll(prevSolution.getCoords());
					for (int j = values.size(); j < numberOfVariables(); j++) {
						values.add(field.zero());
					}
					if (evaluate(p, values).equals(field.zero())) {
						partialSolutions.add(prevSolution);
					}
				}
			} else {
				if (variable > 0) {
					map[variable - 1] = -1;
				}
				map[highest - 1] = 0;
				variable = highest;
				prevPartialSolutions.clear();
				prevPartialSolutions.addAll(partialSolutions);
				partialSolutions.clear();
				for (AffinePoint<T> prevSolution : prevPartialSolutions) {
					List<T> values = new ArrayList<>();
					values.addAll(prevSolution.getCoords());
					for (int j = values.size(); j < numberOfVariables(); j++) {
						values.add(null);
					}
					Polynomial<T> evaluated = partiallyEvaluate(p, values);
					Polynomial<T> univar = field.getUnivariatePolynomialRing()
							.normalize(field.getUnivariatePolynomialRing().getEmbedding(evaluated, map));
					List<T> roots = field.roots(univar);
					for (T root : roots) {
						List<T> coords = new ArrayList<>();
						coords.addAll(prevSolution.getCoords());
						coords.add(root);
						partialSolutions.add(new AffinePoint<>(field, coords));
					}
				}
			}
		}
		if (highest != numberOfVariables()) {
			throw new ArithmeticException("Too many degrees of freedom to solve!");
		}
		return partialSolutions;
	}

	private int highestVariable(Polynomial<T> t) {
		int highest = 0;
		for (Monomial m : t.monomials()) {
			if (t.coefficient(m).equals(getRing().zero())) {
				continue;
			}
			for (int i = numberOfVariables(); i >= 1; i--) {
				if (m.exponents()[i - 1] != 0) {
					highest = Math.max(highest, i);
				}
			}
		}
		return highest;
	}
}
