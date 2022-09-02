package fields.polynomials;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import fields.exceptions.InfinityException;
import fields.helper.AbstractAlgebra;
import fields.helper.AbstractElement;
import fields.interfaces.Element;
import fields.interfaces.Ideal;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.polynomials.DifferentialForms.DifferentialForm;
import fields.vectors.ExteriorProduct;
import fields.vectors.ExteriorProduct.WedgeProduct;
import fields.vectors.ExteriorProduct.WedgeProductSum;
import fields.vectors.FreeModule;
import fields.vectors.Vector;

public class DifferentialForms<T extends Element<T>> extends AbstractAlgebra<Polynomial<T>, DifferentialForm<T>> {
	public static class DifferentialForm<T extends Element<T>> extends AbstractElement<DifferentialForm<T>> {
		private WedgeProductSum<Polynomial<T>> asWedgeProduct;

		private DifferentialForm(WedgeProductSum<Polynomial<T>> asWedgeProduct) {
			this.asWedgeProduct = asWedgeProduct;
		}

		@Override
		public String toString() {
			return asWedgeProduct.toString();
		}

		@Override
		public int compareTo(DifferentialForm<T> o) {
			return asWedgeProduct.compareTo(o.asWedgeProduct);
		}
	}

	private PolynomialRing<T> polynomialRing;
	private ExteriorProduct<Polynomial<T>, Vector<Polynomial<T>>> asExteriorProduct;

	DifferentialForms(PolynomialRing<T> polynomialRing) {
		this.polynomialRing = polynomialRing;
		FreeModule<Polynomial<T>> module = new FreeModule<>(polynomialRing, polynomialRing.numberOfVariables());
		String[] names = new String[polynomialRing.numberOfVariables()];
		for (int i = 0; i < polynomialRing.numberOfVariables(); i++) {
			names[i] = "d" + polynomialRing.getVariableNames()[i];
		}
		this.asExteriorProduct = new ExteriorProduct<>(module, names);
	}

	@Override
	public PolynomialRing<T> getRing() {
		return polynomialRing;
	}

	private DifferentialForm<T> fromWedgeProductSum(WedgeProductSum<Polynomial<T>> t) {
		return new DifferentialForm<>(t);
	}

	private WedgeProductSum<Polynomial<T>> asWedgeProductSum(DifferentialForm<T> t) {
		return t.asWedgeProduct;
	}

	@Override
	public DifferentialForm<T> zero() {
		return fromWedgeProductSum(asExteriorProduct.zero());
	}

	@Override
	public DifferentialForm<T> one() {
		return fromWedgeProductSum(asExteriorProduct.one());
	}

	@Override
	public DifferentialForm<T> add(DifferentialForm<T> s1, DifferentialForm<T> s2) {
		return fromWedgeProductSum(asExteriorProduct.add(asWedgeProductSum(s1), asWedgeProductSum(s2)));
	}

	@Override
	public DifferentialForm<T> negative(DifferentialForm<T> s) {
		return fromWedgeProductSum(asExteriorProduct.negative(asWedgeProductSum(s)));
	}

	@Override
	public DifferentialForm<T> multiply(DifferentialForm<T> s1, DifferentialForm<T> s2) {
		return fromWedgeProductSum(asExteriorProduct.multiply(asWedgeProductSum(s1), asWedgeProductSum(s2)));
	}

	@Override
	public boolean isFree() {
		return true;
	}

	@Override
	public PolynomialIdeal<T> annihilator() {
		return polynomialRing.getZeroIdeal();
	}

	@Override
	public List<Vector<Polynomial<T>>> getSyzygies() {
		return Collections.emptyList();
	}

	private List<WedgeProductSum<Polynomial<T>>> asWedgeProductSumList(List<DifferentialForm<T>> s) {
		List<WedgeProductSum<Polynomial<T>>> result = new ArrayList<>();
		for (DifferentialForm<T> t : s) {
			result.add(asWedgeProductSum(t));
		}
		return result;
	}

	@Override
	public boolean isLinearIndependent(List<DifferentialForm<T>> s) {
		return asExteriorProduct.isLinearIndependent(asWedgeProductSumList(s));
	}

	@Override
	public boolean isGeneratingModule(List<DifferentialForm<T>> s) {
		return asExteriorProduct.isGeneratingModule(asWedgeProductSumList(s));
	}

	@Override
	public List<Vector<Polynomial<T>>> nonTrivialCombinations(List<DifferentialForm<T>> s) {
		return asExteriorProduct.nonTrivialCombinations(asWedgeProductSumList(s));
	}

	private List<DifferentialForm<T>> fromWedgeProductSumList(List<WedgeProductSum<Polynomial<T>>> s) {
		List<DifferentialForm<T>> result = new ArrayList<>();
		for (WedgeProductSum<Polynomial<T>> t : s) {
			result.add(fromWedgeProductSum(t));
		}
		return result;
	}

	@Override
	public List<DifferentialForm<T>> getModuleGenerators() {
		return fromWedgeProductSumList(asExteriorProduct.getModuleGenerators());
	}

	public List<DifferentialForm<T>> getGradedModuleGenerators(int degree) {
		return fromWedgeProductSumList(asExteriorProduct.getGradedModuleGenerators(degree));
	}

	public FreeModule<Polynomial<T>> asFreeModule() {
		return asExteriorProduct.asFreeModule();
	}

	public FreeModule<Polynomial<T>> asGradedFreeModule(int degree) {
		return asExteriorProduct.asGradedFreeModule(degree);
	}

	@Override
	public Vector<Polynomial<T>> asVector(DifferentialForm<T> s) {
		return asExteriorProduct.asVector(asWedgeProductSum(s));
	}

	public Vector<Polynomial<T>> asGradedVector(DifferentialForm<T> s, int degree) {
		return asExteriorProduct.asGradedVector(asWedgeProductSum(s), degree);
	}

	public DifferentialForm<T> getEmbedding(Polynomial<T> t) {
		return fromWedgeProductSum(asExteriorProduct.getEmbedding(t));
	}

	public DifferentialForm<T> derivative(Polynomial<T> t) {
		return derivative(getEmbedding(t));
	}

	public DifferentialForm<T> derivative(DifferentialForm<T> t) {
		DifferentialForm<T> result = zero();
		WedgeProductSum<Polynomial<T>> asWedgeProductSum = asWedgeProductSum(t);
		for (WedgeProduct primitive : asWedgeProductSum.getWedgeProducts()) {
			Polynomial<T> coefficient = asWedgeProductSum.coefficient(primitive);
			for (int i = 0; i < polynomialRing.numberOfVariables(); i++) {
				if (primitive.getVariables().contains(i + 1)) {
					continue;
				}
				WedgeProductSum<Polynomial<T>> derivative = asExteriorProduct.getEmbedding(
						polynomialRing.derivative(coefficient, i + 1),
						asExteriorProduct.getWedgeProduct(Collections.singleton(i + 1)));
				WedgeProductSum<Polynomial<T>> wedge = asExteriorProduct.getEmbedding(polynomialRing.one(), primitive);
				result = add(fromWedgeProductSum(asExteriorProduct.multiply(derivative, wedge)), result);
			}
		}
		return result;
	}

	@Override
	public Exactness exactness() {
		return polynomialRing.exactness();
	}

	@Override
	public DifferentialForm<T> getRandomElement() {
		return fromWedgeProductSum(asExteriorProduct.getRandomElement());
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
	public Iterator<DifferentialForm<T>> iterator() {
		return new Iterator<>() {
			private Iterator<WedgeProductSum<Polynomial<T>>> it = asExteriorProduct.iterator();

			@Override
			public boolean hasNext() {
				return it.hasNext();
			}

			@Override
			public DifferentialForm<T> next() {
				return fromWedgeProductSum(it.next());
			}
		};
	}

	@Override
	public boolean isGeneratingAlgebra(List<DifferentialForm<T>> s) {
		return asExteriorProduct.isGeneratingAlgebra(asWedgeProductSumList(s));
	}

	@Override
	public List<DifferentialForm<T>> getAlgebraGenerators() {
		return fromWedgeProductSumList(asExteriorProduct.getAlgebraGenerators());
	}

	@Override
	public BigInteger characteristic() {
		return asExteriorProduct.characteristic();
	}

	@Override
	public boolean isUnit(DifferentialForm<T> t) {
		return asExteriorProduct.isUnit(asWedgeProductSum(t));
	}

	@Override
	public DifferentialForm<T> inverse(DifferentialForm<T> t) {
		return fromWedgeProductSum(asExteriorProduct.inverse(asWedgeProductSum(t)));
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
		return asExteriorProduct.isIrreducible();
	}

	@Override
	public boolean isZeroDivisor(DifferentialForm<T> t) {
		return asExteriorProduct.isZeroDivisor(asWedgeProductSum(t));
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
	public FactorizationResult<DifferentialForm<T>, DifferentialForm<T>> uniqueFactorization(DifferentialForm<T> t) {
		throw new ArithmeticException("No!");
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
	public boolean isDivisible(DifferentialForm<T> dividend, DifferentialForm<T> divisor) {
		return asExteriorProduct.isDivisible(asWedgeProductSum(dividend), asWedgeProductSum(divisor));
	}

	@Override
	public QuotientAndRemainderResult<DifferentialForm<T>> quotientAndRemainder(DifferentialForm<T> dividend,
			DifferentialForm<T> divisor) {
		QuotientAndRemainderResult<WedgeProductSum<Polynomial<T>>> qr = asExteriorProduct
				.quotientAndRemainder(asWedgeProductSum(dividend), asWedgeProductSum(divisor));
		return new QuotientAndRemainderResult<>(fromWedgeProductSum(qr.getQuotient()),
				fromWedgeProductSum(qr.getRemainder()));
	}

	@Override
	public BigInteger euclidMeasure(DifferentialForm<T> t) {
		throw new ArithmeticException("No!");
	}

	@Override
	public DifferentialForm<T> projectToUnit(DifferentialForm<T> t) {
		return fromWedgeProductSum(asExteriorProduct.projectToUnit(asWedgeProductSum(t)));
	}

	@Override
	public Iterable<DifferentialForm<T>> getUnits() {
		return new Iterable<>() {

			@Override
			public Iterator<DifferentialForm<T>> iterator() {
				Iterator<WedgeProductSum<Polynomial<T>>> it = asExteriorProduct.iterator();

				return new Iterator<>() {

					@Override
					public boolean hasNext() {
						return it.hasNext();
					}

					@Override
					public DifferentialForm<T> next() {
						return fromWedgeProductSum(it.next());
					}
				};
			}
		};
	}

	@Override
	public int krullDimension() {
		return asExteriorProduct.krullDimension();
	}

	@Override
	public IdealResult<DifferentialForm<T>, ?> getIdealWithTransforms(List<DifferentialForm<T>> generators) {
		throw new ArithmeticException("no!");
	}

	@Override
	public Ideal<DifferentialForm<T>> intersect(Ideal<DifferentialForm<T>> t1, Ideal<DifferentialForm<T>> t2) {
		throw new ArithmeticException("no!");
	}

	@Override
	public Ideal<DifferentialForm<T>> radical(Ideal<DifferentialForm<T>> t) {
		throw new ArithmeticException("no!");
	}

}
