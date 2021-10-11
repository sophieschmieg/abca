package varieties.affine;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import fields.interfaces.Element;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.polynomials.CoordinateRing;
import fields.polynomials.Monomial;
import fields.polynomials.PolynomialIdeal;
import fields.polynomials.CoordinateRing.CoordinateIdeal;
import fields.polynomials.CoordinateRing.CoordinateRingElement;
import varieties.Morphism;
import varieties.SpectrumOfField;

public class AffineMorphism<T extends Element<T>> implements Morphism<T, AffinePoint<T>, AffinePoint<T>> {
	private AffineScheme<T> domain;
	private AffineScheme<T> range;
	private List<Polynomial<T>> polynomials;

	public AffineMorphism(AffineScheme<T> domain, AffineScheme<T> range, List<CoordinateRingElement<T>> polynomials) {
		this.domain = domain;
		this.range = range;
		this.polynomials = new ArrayList<>();
		for (CoordinateRingElement<T> p : polynomials) {
			this.polynomials.add(p.getElement());
		}
		PolynomialRing<T> domainRing = domain.getCoordinateRing().getPolynomialRing();
		for (Polynomial<T> rangeGenerator : range.getCoordinateRing().getIdeal().generators()) {
			Polynomial<T> substitued = domainRing.substitute(rangeGenerator, this.polynomials);
			if (!domain.getCoordinateRing().getIdeal().contains(substitued)) {
				throw new ArithmeticException("Affine morphism ill defined");
			}
		}
	}

	public static <T extends Element<T>> AffineMorphism<T> fromPolynomials(AffineScheme<T> domain,
			AffineScheme<T> range, List<Polynomial<T>> polynomials) {
		List<CoordinateRingElement<T>> result = new ArrayList<>();
		for (Polynomial<T> p : polynomials) {
			result.add(domain.getCoordinateRing().getEmbedding(p));
		}
		return new AffineMorphism<>(domain, range, result);
	}

	public static <T extends Element<T>> AffineMorphism<T> composition(AffineMorphism<T> first,
			AffineMorphism<T> second) {
		List<Polynomial<T>> polynomials = new ArrayList<>();
		PolynomialRing<T> domainRing = first.domain.getCoordinateRing().getPolynomialRing();
		for (Polynomial<T> secondPolynomial : second.polynomials) {
			polynomials.add(domainRing.substitute(secondPolynomial, first.polynomials));
		}
		return fromPolynomials(first.domain, second.range, polynomials);
	}

	public static class AffineFiberedProduct<T extends Element<T>> {
		private AffineScheme<T> product;
		private AffineScheme<T> factor1;
		private AffineScheme<T> factor2;
		private AffineScheme<T> base;
		private AffineMorphism<T> projection1;
		private AffineMorphism<T> projection2;
		private AffineMorphism<T> structural1;
		private AffineMorphism<T> structural2;

		public AffineScheme<T> getProduct() {
			return product;
		}

		public AffineScheme<T> getFactor1() {
			return factor1;
		}

		public AffineScheme<T> getFactor2() {
			return factor2;
		}

		public AffineScheme<T> getBase() {
			return base;
		}

		public AffineMorphism<T> getProjection1() {
			return projection1;
		}

		public AffineMorphism<T> getProjection2() {
			return projection2;
		}

		public AffineMorphism<T> getStructural1() {
			return structural1;
		}

		public AffineMorphism<T> getStructural2() {
			return structural2;
		}
	}

	public static <T extends Element<T>> AffineFiberedProduct<T> fiberedProduct(AffineMorphism<T> structural1,
			AffineMorphism<T> structural2) {
		if (!structural1.getRange().equals(structural2.getRange())) {
			throw new ArithmeticException("Not defined over the same base");
		}

		AffineFiberedProduct<T> result = new AffineFiberedProduct<>();
		result.structural1 = structural1;
		result.structural2 = structural2;
		result.factor1 = structural1.getDomain();
		result.factor2 = structural2.getDomain();
		result.base = structural1.getRange();
		CoordinateRing<T> coordinateRing1 = result.factor1.getCoordinateRing();
		PolynomialRing<T> polynomialRing1 = coordinateRing1.getPolynomialRing();
		CoordinateRing<T> coordinateRing2 = result.factor2.getCoordinateRing();
		PolynomialRing<T> polynomialRing2 = coordinateRing2.getPolynomialRing();

		CoordinateRing<T> tensorProduct = CoordinateRing.tensorProduct(coordinateRing1, coordinateRing2);
		PolynomialRing<T> polynomialRing = tensorProduct.getPolynomialRing();
		int[] map1 = new int[polynomialRing1.numberOfVariables()];
		for (int i = 0; i < map1.length; i++) {
			map1[i] = i;
		}
		int[] map2 = new int[polynomialRing2.numberOfVariables()];
		for (int i = 0; i < map2.length; i++) {
			map2[i] = i + map1.length;
		}
		List<CoordinateRingElement<T>> additionalEquations = new ArrayList<>();
		for (int i = 0; i < result.base.getCoordinateRing().getPolynomialRing().numberOfVariables(); i++) {
			Polynomial<T> x1 = polynomialRing.getEmbedding(structural1.polynomials.get(i), map1);
			Polynomial<T> x2 = polynomialRing.getEmbedding(structural2.polynomials.get(i), map2);
			Polynomial<T> equation = polynomialRing.subtract(x1, x2);
			additionalEquations.add(tensorProduct.getEmbedding(equation));
		}
		CoordinateIdeal<T> ideal = tensorProduct.getIdeal(additionalEquations);
		CoordinateRing<T> pushout = ideal.divideOut();
		result.product = new AffineScheme<>(result.factor1.getField(), pushout);
		List<CoordinateRingElement<T>> projection1List = new ArrayList<>();
		for (int i = 0; i < polynomialRing1.numberOfVariables(); i++) {
			projection1List.add(tensorProduct.getVar(i + 1));
		}
		List<CoordinateRingElement<T>> projection2List = new ArrayList<>();
		for (int i = 0; i < polynomialRing2.numberOfVariables(); i++) {
			projection2List.add(tensorProduct.getVar(i + 1 + polynomialRing1.numberOfVariables()));
		}
		result.projection1 = new AffineMorphism<>(result.product, result.factor1, projection1List);
		result.projection2 = new AffineMorphism<>(result.product, result.factor2, projection2List);
		return result;
	}

	@Override
	public String toString() {
		StringBuilder build = new StringBuilder();
		build.append("(");
		boolean first = true;
		for (Polynomial<T> polynomial : polynomials) {
			if (first) {
				first = false;
			} else {
				build.append(", ");
			}
			build.append(polynomial);
		}
		build.append(")");
		return build.toString();
	}

	@Override
	public AffinePoint<T> evaluate(AffinePoint<T> t) {
		List<T> list = new ArrayList<>();
		for (Polynomial<T> polynomial : polynomials) {
			list.add(domain.getCoordinateRing().getPolynomialRing().evaluate(polynomial, t));
		}
		return new AffinePoint<>(domain.getField(), list);
	}

	public AffineScheme<T> getDomain() {
		return domain;
	}

	public AffineScheme<T> getRange() {
		return range;
	}

	public List<Polynomial<T>> getPolynomials() {
		return polynomials;
	}

	public boolean isInjective() {
		CoordinateRing<T> domainRing = domain.getCoordinateRing();
		List<CoordinateRingElement<T>> cres = new ArrayList<>();
		for (Polynomial<T> polynomial : polynomials) {
			cres.add(domainRing.getEmbedding(polynomial));
		}
		CoordinateIdeal<T> ideal = domainRing.getIdeal(cres);
		for (CoordinateRingElement<T> generator : ideal.generators()) {
			if (generator.getElement().degree() > 1) {
				return false;
			}
		}
		return ideal.divideOut().krullDimension() <= 0;
	}

	public boolean hasDenseImage() {
		AffineMorphism<T> image = image();
		if (image.getDomain().dimension() != range.dimension()) {
			return false;
		}
		return true; // maybe
	}

	@Override
	public RestrictionResult<T> restrict(AffinePoint<T> p) {
		return new RestrictionResult<>(0, 0, domain.identityMorphism(), range.identityMorphism(), this);
	}

	public AffineMorphism<T> image() {
		CoordinateRing<T> domainCoordinateRing = domain.getCoordinateRing();
		PolynomialRing<T> domainPolynomialRing = domainCoordinateRing.getPolynomialRing();
		CoordinateRing<T> rangeCoordinateRing = range.getCoordinateRing();
		PolynomialRing<T> rangePolynomialRing = rangeCoordinateRing.getPolynomialRing();
		Comparator<Monomial> eliminationOrder = new Monomial.EliminationOrder(Monomial.GREVLEX, Monomial.GREVLEX,
				domainPolynomialRing.numberOfVariables());
		CoordinateRing<T> tensorProduct = CoordinateRing.tensorProduct(domainCoordinateRing, rangeCoordinateRing,
				eliminationOrder);
		PolynomialRing<T> productPolynomialRing = tensorProduct.getPolynomialRing();
		List<CoordinateRingElement<T>> graphConstraints = new ArrayList<>();
		for (int i = 0; i < rangePolynomialRing.numberOfVariables(); i++) {
			Polynomial<T> polynomial = polynomials.get(i);
			Polynomial<T> inProduct = productPolynomialRing.getEmbedding(polynomial);
			Polynomial<T> constraint = productPolynomialRing.subtract(inProduct,
					productPolynomialRing.getVar(i + 1 + domainPolynomialRing.numberOfVariables()));
			graphConstraints.add(tensorProduct.getEmbedding(constraint));
		}
		CoordinateIdeal<T> graphIdeal = tensorProduct.getIdeal(graphConstraints);
		PolynomialIdeal<T> eliminatedIdeal = graphIdeal.asPolynomialIdeal();
		List<CoordinateRingElement<T>> imageGenerators = new ArrayList<>();
		int map[] = new int[productPolynomialRing.numberOfVariables()];
		for (int i = 0; i < domainPolynomialRing.numberOfVariables(); i++) {
			map[i] = -1;
		}
		List<CoordinateRingElement<T>> embedding = new ArrayList<>();
		for (int i = 0; i < rangePolynomialRing.numberOfVariables(); i++) {
			map[i + domainPolynomialRing.numberOfVariables()] = i;
			embedding.add(rangeCoordinateRing.getVar(i + 1));
		}
		generatorLoop: for (Polynomial<T> graphGenerator : eliminatedIdeal.generators()) {
			Monomial leading = graphGenerator.leadingMonomial();
			for (int i = 0; i < domainPolynomialRing.numberOfVariables(); i++) {
				if (leading.exponents()[i] > 0) {
					continue generatorLoop;
				}
			}
			imageGenerators
					.add(rangeCoordinateRing.getEmbedding(rangePolynomialRing.getEmbedding(graphGenerator, map)));
		}
		AffineScheme<T> image = new AffineScheme<>(range.getField(),
				rangeCoordinateRing.getIdeal(imageGenerators).divideOut());
		return new AffineMorphism<>(image, range, embedding);
	}

	/**
	 * Returns the closure of the preimage given by an injection into the range.
	 * 
	 * @return
	 */

	public AffineMorphism<T> preimage(AffineMorphism<T> image) {
		if (!image.isInjective()) {
			throw new ArithmeticException("image not an injection into range!");
		}
		return fiberedProduct(this, image).projection1;
	}

	public AffineMorphism<T> preimage(AffinePoint<T> p) {
		return preimage(range.pointAsMorphism(p).restrict(SpectrumOfField.POINT).getRestrictedMorphism());
	}
}
