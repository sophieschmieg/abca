package fields.polynomials;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.helper.AbstractIdeal;
import fields.interfaces.Element;
import fields.interfaces.Ideal;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.Ring.PrimaryDecompositionResult;
import fields.vectors.FreeModule;
import fields.vectors.Vector;

public class PolynomialIdeal<T extends Element<T>> extends AbstractIdeal<Polynomial<T>> {
	private List<Polynomial<T>> basis;
	private List<Vector<Polynomial<T>>> syzygies;
	private PolynomialRing<T> polynomialRing;
	private CoordinateRing<T> coordinateRing;

	PolynomialIdeal(PolynomialRing<T> polynomialRing, List<Polynomial<T>> generators) {
		super(polynomialRing);
		this.polynomialRing = polynomialRing;
		this.basis = generators;
	}

	@Override
	public PolynomialRing<T> getRing() {
		return polynomialRing;
	}

	public CoordinateRing<T> divideOut() {
		if (coordinateRing == null) {
			coordinateRing = new CoordinateRing<>(polynomialRing, this);
		}
		return coordinateRing;
	}

	public int dimension() {
		return divideOut().krullDimension() - polynomialRing.getRing().krullDimension();
	}

	public int degree() {
		return divideOut().degree();
	}

//	@Override
//	public List<Vector<Polynomial<T>>> nonTrivialCombinations(List<Polynomial<T>> s) {
//		Ring<T> ring = polynomialRing.getRing();
//		IdealResult<Polynomial<T>, PolynomialIdeal<T>> ideal = polynomialRing.getIdealWithTransforms(s);
//		if (ideal.getIdeal().generators().equals(s)) {
//			return ideal.getIdeal().getModuleGeneratorRelations();
//		}
//		List<List<Polynomial<T>>> result = new ArrayList<>();
//		for (int i = 1; i < s.size(); i++) {
//			for (int j = 0; j < i; j++) {
//				Monomial leadingFirst = s.get(j).leadingMonomial();
//				T firstLeadingCoefficient = s.get(j).leadingCoefficient();
//				Monomial leadingSecond = s.get(i).leadingMonomial();
//				T secondLeadingCoefficient = s.get(i).leadingCoefficient();
//				Monomial lcm = leadingFirst.lcm(leadingSecond);
//				T coefficientLcm = ring.lcm(firstLeadingCoefficient, secondLeadingCoefficient);
//				Polynomial<T> firstMultiplier = polynomialRing.getPolynomial(Collections.singletonMap(
//						lcm.divide(leadingFirst), ring.divideChecked(coefficientLcm, firstLeadingCoefficient)));
//				Polynomial<T> secondMultiplier = polynomialRing
//						.getPolynomial(Collections.singletonMap(lcm.divide(leadingSecond),
//								ring.negative(ring.divideChecked(coefficientLcm, secondLeadingCoefficient))));
//				Polynomial<T> buchbergerPolynomial = polynomialRing.add(
//						polynomialRing.multiply(firstMultiplier, s.get(j)),
//						polynomialRing.multiply(secondMultiplier, s.get(i)));
//				List<Polynomial<T>> generated = ideal.getIdeal().generate(buchbergerPolynomial);
//				List<Polynomial<T>> relation = new ArrayList<>();
//				for (int k = 0; k < s.size(); k++) {
//					if (k == j) {
//						relation.add(firstMultiplier);
//					} else if (k == i) {
//						relation.add(secondMultiplier);
//					} else {
//						relation.add(polynomialRing.zero());
//					}
//				}
//				for (int k = 0; k < generated.size(); k++) {
//					Polynomial<T> coefficient = polynomialRing.negative(generated.get(k));
//					for (int l = 0; l < s.size(); l++) {
//						Polynomial<T> originalCoefficient = polynomialRing.multiply(coefficient,
//								ideal.getGeneratorExpressions().get(k).get(l));
//						relation.set(l, polynomialRing.add(relation.get(l), originalCoefficient));
//					}
//				}
//				result.add(relation);
//			}
//		}
//		return result;
//	}

//	@Override
//	public List<List<Polynomial<T>>> getModuleGeneratorRelations() {
//		List<List<Polynomial<T>>> result = new ArrayList<>();
//		for (Vector<Polynomial<T>> syzygy : syzygies()) {
//			result.add(syzygy.asList());
//		}
//		return result;
//	}

	@Override
	public boolean contains(Polynomial<T> polynomial) {
		return polynomialRing.reduce(polynomial, this.basis).equals(polynomialRing.zero());
	}

	@Override
	public Polynomial<T> residue(Polynomial<T> f) {
		return polynomialRing.reduce(f, this.basis);
	}

	public PolynomialIdeal<T> saturate(Polynomial<T> by) {
		PolynomialRing<T> ring = polynomialRing.addVariableWithElimination(1);
		Polynomial<T> t = ring.getVar(1);
		Polynomial<T> f = ring.getEmbeddingWithElimination(by, 1);
		Polynomial<T> inv = ring.subtract(ring.one(), ring.multiply(f, t));
		List<Polynomial<T>> generators = new ArrayList<Polynomial<T>>();
		generators.add(inv);
		for (Polynomial<T> p : this.basis) {
			generators.add(ring.getEmbeddingWithElimination(p, 1));
		}
		Ideal<Polynomial<T>> saturationIdeal = ring.getIdeal(generators);
		List<Polynomial<T>> saturationGenerators = new ArrayList<Polynomial<T>>();
		for (Polynomial<T> b : saturationIdeal.generators()) {
			if (b.leadingMonomial().exponents()[0] == 0)
				saturationGenerators.add(polynomialRing.getEmbeddingWithElimination(b, -1));
		}
		return new PolynomialIdeal<T>(polynomialRing, saturationGenerators);
	}

	public PolynomialIdeal<T> saturate(PolynomialIdeal<T> by) {
		PolynomialIdeal<T> result = this;
		for (Polynomial<T> generator : by.generators()) {
			result = result.saturate(generator);
		}
		return result;
	}

	public String toString() {
		return this.basis.toString();
	}

	public boolean isFinite() {
		return false;
	}

	@Override
	public BigInteger getNumberOfElements() {
		throw new InfinityException();
	}

	@Override
	public List<Polynomial<T>> generators() {
		return Collections.unmodifiableList(this.basis);
	}

	@Override
	public List<Polynomial<T>> generate(Polynomial<T> t) {
		return polynomialRing.generalQuotientAndRemainder(t, basis).getQuotients();
	}

	public Ideal<T> intersectToRing() {
		List<T> ringGenerators = new ArrayList<>();
		for (Polynomial<T> generator : generators()) {
			if (generator.degree() == 0) {
				ringGenerators.add(generator.leadingCoefficient());
			}
		}
		return polynomialRing.getRing().getIdeal(ringGenerators);
	}

	public boolean isHomogenous() {
		for (Polynomial<T> generator : generators()) {
			if (!polynomialRing.isHomogeneous(generator)) {
				return false;
			}
		}
		return true;
	}

	@Override
	public boolean isPrimary() {
		return polynomialRing.primaryDecomposition(this).getPrimaryIdeals().size() == 1;
	}

	@Override
	public boolean isPrime() {
		PrimaryDecompositionResult<Polynomial<T>, PolynomialIdeal<T>> decomposition = polynomialRing
				.primaryDecomposition(this);
		if (decomposition.getPrimaryIdeals().size() != 1) {
			return false;
		}
		return decomposition.getPrimaryIdeals().get(0).equals(decomposition.getRadicals().get(0));
	}

	@Override
	public boolean isMaximal() {
		return isPrime() && divideOut().krullDimension() == polynomialRing.getRing().krullDimension()
				&& intersectToRing().isMaximal();
	}

	@Override
	public boolean isRadical() {
		PrimaryDecompositionResult<Polynomial<T>, PolynomialIdeal<T>> decomposition = polynomialRing
				.primaryDecomposition(this);
		for (int i = 0; i < decomposition.getPrimaryIdeals().size(); i++) {
			if (!decomposition.getPrimaryIdeals().get(i).equals(decomposition.getRadicals().get(i))) {
				return false;
			}
		}
		return true;
	}

	public List<PolynomialIdeal<T>> minimalPrimeIdealsOver() {
		return polynomialRing.primaryDecomposition(this).getRadicals();
	}

	public boolean isMaximalOverAlgebraicClosure() {
		if (polynomialRing.numberOfVariables() != polynomialRing.krullDimension()) {
			throw new ArithmeticException("Only defined over field!");
		}
		if (polynomialRing.numberOfVariables() != generators().size()) {
			return false;
		}
		for (Polynomial<T> generator : generators()) {
			if (generator.degree() != 1) {
				return false;
			}
		}
		return true;
	}

	public boolean isCompleteIntersection() {
		return dimension() == polynomialRing.numberOfVariables() - generators().size();
	}

	public List<Vector<Polynomial<T>>> syzygies() {
		if (basis.size() == 1 && basis.get(0).equals(polynomialRing.zero())) {
			return Collections.singletonList(new Vector<>(polynomialRing.one()));
		}
		if (syzygies == null) {
			FreeModule<Polynomial<T>> module = new FreeModule<>(polynomialRing, basis.size());
			syzygies = new ArrayList<>();
			for (int i = 0; i < basis.size(); i++) {
				for (int j = 0; j < i; j++) {
					Monomial mi = basis.get(i).leadingMonomial();
					Monomial mj = basis.get(j).leadingMonomial();
					Monomial mlcm = mi.lcm(mj);
					Monomial lmi = mlcm.divide(mi);
					Monomial lmj = mlcm.divide(mj);
					T ci = basis.get(i).leadingCoefficient();
					T cj = basis.get(j).leadingCoefficient();
					T clcm = polynomialRing.getRing().lcm(ci, cj);
					T lci = polynomialRing.getRing().divideChecked(clcm, ci);
					T lcj = polynomialRing.getRing().negative(polynomialRing.getRing().divideChecked(clcm, cj));
					Polynomial<T> cofactorI = polynomialRing.getEmbedding(lci, lmi.exponents());
					Polynomial<T> cofactorJ = polynomialRing.getEmbedding(lcj, lmj.exponents());
					Polynomial<T> combined = polynomialRing.add(polynomialRing.multiply(cofactorI, basis.get(i)),
							polynomialRing.multiply(cofactorJ, basis.get(j)));
					List<Polynomial<T>> generate = generate(combined);
					List<Polynomial<T>> syzygy = new ArrayList<>();
					syzygy.addAll(generate);
					syzygy.set(i, polynomialRing.subtract(syzygy.get(i), cofactorI));
					syzygy.set(j, polynomialRing.subtract(syzygy.get(j), cofactorJ));
					Vector<Polynomial<T>> syzygyVector = new Vector<>(syzygy);
					if (!syzygyVector.equals(module.zero())) {
						syzygies.add(syzygyVector);
					}
				}
			}
		}
		return syzygies;
	}

}
