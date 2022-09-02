package varieties.affine;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.TreeSet;

import fields.helper.TranscendentalFieldExtension;
import fields.helper.TranscendentalFieldExtension.TExt;
import fields.interfaces.Element;
import fields.interfaces.Polynomial;
import fields.interfaces.PolynomialRing;
import fields.interfaces.UnivariatePolynomial;
import fields.interfaces.UnivariatePolynomialRing;
import fields.polynomials.AbstractPolynomialRing;
import fields.polynomials.CoordinateRing;
import fields.polynomials.CoordinateRing.CoordinateIdeal;
import fields.polynomials.CoordinateRing.CoordinateRingElement;
import fields.polynomials.Monomial;
import fields.polynomials.PolynomialIdeal;
import varieties.AbstractMorphism;
import varieties.GeneralRationalFunction;
import varieties.SpectrumOfField;
import varieties.affine.AffineScheme.AffineSimplification;

public class AffineMorphism<T extends Element<T>> extends AbstractMorphism<T, AffinePoint<T>, AffinePoint<T>> {
	private AffineScheme<T> domain;
	private AffineScheme<T> range;
	private List<Polynomial<T>> polynomials;
	private AffineImmersionClassification<T> immersionClassification;

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
				throw new ArithmeticException("Affine morphism ill defined: domain: " + domain + " range: " + range
						+ " polynomials: " + polynomials + " substituted: " + substitued);
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

	public static <T extends Element<T>> AffineMorphism<T> getClosedImmersion(AffineScheme<T> range,
			CoordinateIdeal<T> ideal) {
		AffineScheme<T> domain = new AffineScheme<>(range.getField(), ideal.divideOut());
		List<CoordinateRingElement<T>> morphism = new ArrayList<>();
		for (int i = 0; i < range.getCoordinateRing().getPolynomialRing().numberOfVariables(); i++) {
			morphism.add(domain.getCoordinateRing().getVar(i + 1));
		}
		return new AffineMorphism<>(domain, range, morphism).simplifyDomain();
	}

	public static class OpenImmersionResult<T extends Element<T>> {
		private AffineMorphism<T> openImmersion;
		private CoordinateRingElement<T> polynomial;
		private CoordinateRingElement<T> inverse;

		private OpenImmersionResult(AffineMorphism<T> openImmersion, CoordinateRingElement<T> polynomial,
				CoordinateRingElement<T> inverse) {
			this.openImmersion = openImmersion;
			this.polynomial = polynomial;
			this.inverse = inverse;
		}

		public AffineMorphism<T> getOpenImmersion() {
			return openImmersion;
		}

		public CoordinateRingElement<T> getPolynomial() {
			return polynomial;
		}

		public CoordinateRingElement<T> getInverse() {
			return inverse;
		}
	}

	public static <T extends Element<T>> OpenImmersionResult<T> getOpenImmersion(AffineScheme<T> range,
			Polynomial<T> polynomial) {
		polynomial = range.getCoordinateRing().getEmbedding(polynomial).getElement();
		polynomial = range.getCoordinateRing().getPolynomialRing().radical(polynomial);
		PolynomialRing<T> domainPolynomialRing = AbstractPolynomialRing.getPolynomialRing(range.getField(),
				range.getCoordinateRing().getPolynomialRing().numberOfVariables() + 1, Monomial.GREVLEX);
		List<Polynomial<T>> generators = new ArrayList<>();
		for (Polynomial<T> generator : range.getCoordinateRing().getIdeal().generators()) {
			generators.add(domainPolynomialRing.getEmbedding(generator));
		}
		Polynomial<T> generator = domainPolynomialRing.subtract(
				domainPolynomialRing.multiply(domainPolynomialRing.getEmbedding(polynomial),
						domainPolynomialRing.getVar(domainPolynomialRing.numberOfVariables())),
				domainPolynomialRing.one());
		generators.add(generator);
		AffineScheme<T> domain = new AffineScheme<>(range.getField(),
				domainPolynomialRing.getIdeal(generators).divideOut());
		List<CoordinateRingElement<T>> morphism = new ArrayList<>();
		for (int i = 0; i < range.getCoordinateRing().getPolynomialRing().numberOfVariables(); i++) {
			morphism.add(domain.getCoordinateRing().getVar(i + 1));
		}
		AffineMorphism<T> immersion = new AffineMorphism<>(domain, range, morphism);
		AffineSimplification<T> simplificiation = immersion.getDomain().simplify();
		return new OpenImmersionResult<>(AffineMorphism.composition(simplificiation.getInverse(), immersion),
				range.getCoordinateRing().getEmbedding(polynomial), simplificiation.getInverse()
						.inducedMap(domain.getCoordinateRing().getVar(domainPolynomialRing.numberOfVariables())));
	}

	public static <T extends Element<T>> OpenImmersionResult<T> getOpenImmersion(AffineScheme<T> range,
			List<Polynomial<T>> polynomials) {
		PolynomialRing<T> rangePolynomialRing = range.getCoordinateRing().getPolynomialRing();
		Polynomial<T> polynomial = rangePolynomialRing.one();
		for (Polynomial<T> constraint : polynomials) {
			polynomial = rangePolynomialRing.lcm(constraint, polynomial);
		}
		return getOpenImmersion(range, polynomial);
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
		private AffineSimplification<T> simplification;
		private AffineScheme<T> simplifiedProduct;
		private AffineMorphism<T> simplifiedProjection1;
		private AffineMorphism<T> simplifiedProjection2;

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

		public AffineMorphism<T> universalProperty(AffineMorphism<T> projection1, AffineMorphism<T> projection2) {
			if (!AffineMorphism.composition(projection1, structural1)
					.equals(AffineMorphism.composition(projection2, structural2))) {
				throw new ArithmeticException("Diagram does not commute!");
			}
			List<Polynomial<T>> polynomials = new ArrayList<>();
			polynomials.addAll(projection1.getPolynomials());
			polynomials.addAll(projection2.getPolynomials());
			return AffineMorphism.fromPolynomials(projection1.getDomain(), product, polynomials);
		}

		public AffineScheme<T> getSimplifiedProduct() {
			if (simplifiedProduct == null) {
				simplifiedProduct = simplification().getSimplification().getRange();
			}
			return simplifiedProduct;
		}

		public AffineMorphism<T> getSimplifiedProjection1() {
			if (simplifiedProjection1 == null) {
				simplifiedProjection1 = AffineMorphism.composition(simplification().getInverse(), projection1);
			}
			return simplifiedProjection1;
		}

		public AffineMorphism<T> getSimplifiedProjection2() {
			if (simplifiedProjection2 == null) {
				simplifiedProjection2 = AffineMorphism.composition(simplification().getInverse(), projection2);
			}
			return simplifiedProjection2;
		}

		public AffineMorphism<T> simplifiedUniversalProperty(AffineMorphism<T> projection1,
				AffineMorphism<T> projection2) {
			return AffineMorphism.composition(universalProperty(projection1, projection2),
					simplification().getSimplification());
		}

		private AffineSimplification<T> simplification() {
			if (simplification == null) {
				simplification = product.simplify();
			}
			return simplification;
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

	public CoordinateRingElement<T> inducedMap(CoordinateRingElement<T> rangeElement) {
		Polynomial<T> domainPolynomial = domain.getCoordinateRing().getPolynomialRing()
				.substitute(rangeElement.getElement(), polynomials);
		return domain.getCoordinateRing().getEmbedding(domainPolynomial);
	}

	public Polynomial<T> inducedMap(Polynomial<T> rangeElement) {
		return inducedMap(range.getCoordinateRing().getEmbedding(rangeElement)).getElement();
	}
//	
//	public static class IntersectionResult<T extends Element<T>> {
//		private AffineScheme<T> intersection;
//		private AffineMorphism<T> firstEmbedding;
//		private AffineMorphism<T> secondEmbedding;
//
//		private IntersectionResult(AffineScheme<T> intersection, AffineMorphism<T> firstEmbedding,
//				AffineMorphism<T> secondEmbedding) {
//			this.intersection = intersection;
//			this.firstEmbedding = firstEmbedding;
//			this.secondEmbedding = secondEmbedding;
//		}
//
//		public AffineScheme<T> getIntersection() {
//			return intersection;
//		}
//
//		public AffineMorphism<T> getFirstEmbedding() {
//			return firstEmbedding;
//		}
//
//		public AffineMorphism<T> getSecondEmbedding() {
//			return secondEmbedding;
//		}
//	}
//
//	public static <T extends Element<T>> IntersectionResult<T> intersectClosedImmersions(AffineMorphism<T> closed1,
//			AffineMorphism<T> closed2) {
//		if (!closed1.isClosedImmersion() || !closed2.isClosedImmersion()) {
//			throw new ArithmeticException("Not closed immersions!");
//		}
//		if (!closed1.getRange().equals(closed2.getRange())) {
//			throw new ArithmeticException("different ranges!");
//		}
//		
//		CoordinateRing<T> rangeCoordinateRing = closed1.getRange().getCoordinateRing();
//		PolynomialRing<T> rangePolynomialRing = rangeCoordinateRing.getPolynomialRing();
//		CoordinateIdeal<T> ideal = (CoordinateIdeal<T>)rangeCoordinateRing.add(closed1.closedImmersionAsIdeal(),
//				closed2.closedImmersionAsIdeal());
//	AffineMorphism<T> closed =	getClosedImmersion(closed1.getRange(), ideal);
//	CoordinateRing<T> domain1CoordinateRing = closed1.getDomain().getCoordinateRing();
//	PolynomialRing<T> domain1PolynomialRing = domain1CoordinateRing.getPolynomialRing();
//	List<CoordinateRingElement<T>> embedding1 = new ArrayList<>();
//	for (int i = 0; i < domain1PolynomialRing.numberOfVariables(); i++) {
//		TExt<T> inRange = closed1.birationalInverseInducedMap(domain1CoordinateRing.getVar(i + 1));
//		CoordinateRingElement<T> numerator = closed
//				.inducedMap(rangeCoordinateRing.getEmbedding(inRange.getNumerator()));
//		CoordinateRingElement<T> denominator = closed
//				.inducedMap(rangeCoordinateRing.getEmbedding(inRange.getDenominator()));
//		embedding1.add(closed.getDomain().getCoordinateRing().divide(numerator, denominator));
//	}
//	CoordinateRing<T> domain2CoordinateRing = closed2.getDomain().getCoordinateRing();
//	PolynomialRing<T> domain2PolynomialRing = domain2CoordinateRing.getPolynomialRing();
//	List<CoordinateRingElement<T>> embedding2 = new ArrayList<>();
//	for (int i = 0; i < domain2PolynomialRing.numberOfVariables(); i++) {
//		TExt<T> inRange = closed2.birationalInverseInducedMap(domain2CoordinateRing.getVar(i + 1));
//		CoordinateRingElement<T> numerator = open
//				.inducedMap(rangeCoordinateRing.getEmbedding(inRange.getNumerator()));
//		CoordinateRingElement<T> denominator = open
//				.inducedMap(rangeCoordinateRing.getEmbedding(inRange.getDenominator()));
//		embedding2.add(open.getDomain().getCoordinateRing().divide(numerator, denominator));
//	}
//	return new IntersectionResult<>(open.getDomain(),
//			new AffineMorphism<>(open.getDomain(), open1.getDomain(), embedding1),
//			new AffineMorphism<>(open.getDomain(), open2.getDomain(), embedding2));
//	}
//
//	public static <T extends Element<T>> IntersectionResult<T> intersectOpenImmersions(AffineMorphism<T> open1,
//			AffineMorphism<T> open2) {
//		if (!open1.isOpenImmersion() || !open2.isOpenImmersion()) {
//			throw new ArithmeticException("Not open immersions!");
//		}
//		if (!open1.getRange().equals(open2.getRange())) {
//			throw new ArithmeticException("different ranges!");
//		}
//		CoordinateRing<T> rangeCoordinateRing = open1.getRange().getCoordinateRing();
//		CoordinateRingElement<T> combinedObstruction = rangeCoordinateRing.multiply(open1.getInversionObstruction(),
//				open2.getInversionObstruction());
//		AffineMorphism<T> open = getOpenImmersion(open1.getRange(), combinedObstruction.getElement())
//				.getOpenImmersion();
//		CoordinateRing<T> domain1CoordinateRing = open1.getDomain().getCoordinateRing();
//		PolynomialRing<T> domain1PolynomialRing = domain1CoordinateRing.getPolynomialRing();
//		List<CoordinateRingElement<T>> embedding1 = new ArrayList<>();
//		for (int i = 0; i < domain1PolynomialRing.numberOfVariables(); i++) {
//			TExt<T> inRange = open1.birationalInverseInducedMap(domain1CoordinateRing.getVar(i + 1));
//			CoordinateRingElement<T> numerator = open
//					.inducedMap(rangeCoordinateRing.getEmbedding(inRange.getNumerator()));
//			CoordinateRingElement<T> denominator = open
//					.inducedMap(rangeCoordinateRing.getEmbedding(inRange.getDenominator()));
//			embedding1.add(open.getDomain().getCoordinateRing().divide(numerator, denominator));
//		}
//		CoordinateRing<T> domain2CoordinateRing = open2.getDomain().getCoordinateRing();
//		PolynomialRing<T> domain2PolynomialRing = domain2CoordinateRing.getPolynomialRing();
//		List<CoordinateRingElement<T>> embedding2 = new ArrayList<>();
//		for (int i = 0; i < domain2PolynomialRing.numberOfVariables(); i++) {
//			TExt<T> inRange = open2.birationalInverseInducedMap(domain2CoordinateRing.getVar(i + 1));
//			CoordinateRingElement<T> numerator = open
//					.inducedMap(rangeCoordinateRing.getEmbedding(inRange.getNumerator()));
//			CoordinateRingElement<T> denominator = open
//					.inducedMap(rangeCoordinateRing.getEmbedding(inRange.getDenominator()));
//			embedding2.add(open.getDomain().getCoordinateRing().divide(numerator, denominator));
//		}
//		return new IntersectionResult<>(open.getDomain(),
//				new AffineMorphism<>(open.getDomain(), open1.getDomain(), embedding1),
//				new AffineMorphism<>(open.getDomain(), open2.getDomain(), embedding2));
//	}

	@Override
	public String toString() {
		StringBuilder build = new StringBuilder();
		build.append(domain);
		build.append(" --(");
		boolean first = true;
		for (Polynomial<T> polynomial : polynomials) {
			if (first) {
				first = false;
			} else {
				build.append(", ");
			}
			build.append(polynomial);
		}
		build.append(")--> ");
		build.append(range);
		return build.toString();
	}

	@Override
	public boolean equals(Object obj) {
		if (!(obj instanceof AffineMorphism<?>)) {
			return false;
		}
		@SuppressWarnings("unchecked")
		AffineMorphism<T> other = (AffineMorphism<T>) obj;
		if (!domain.equals(other.getDomain())) {
			return false;
		}
		if (!range.equals(other.getRange())) {
			return false;
		}
		return polynomials.equals(other.getPolynomials());
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

	public boolean isFinite() {
		CoordinateRing<T> domainCoordinateRing = domain.getCoordinateRing();
		PolynomialRing<T> domainPolynomialRing = domainCoordinateRing.getPolynomialRing();
		CoordinateRing<T> rangeCoordinateRing = range.getCoordinateRing();
		PolynomialRing<T> rangePolynomialRing = rangeCoordinateRing.getPolynomialRing();
		AffineFiberedProduct<T> graph = graph();
		Comparator<Monomial> eliminationOrder = new Monomial.EliminationOrder(Monomial.REVLEX, Monomial.GREVLEX,
				domainPolynomialRing.numberOfVariables());
		PolynomialRing<T> graphPolynomialRing = AbstractPolynomialRing.getPolynomialRing(domain.getField(),
				domainPolynomialRing.numberOfVariables() + rangePolynomialRing.numberOfVariables(), eliminationOrder);
		PolynomialIdeal<T> graphIdeal = graphPolynomialRing.getEmbedding(graph.product.getCoordinateRing().getIdeal());
		Set<Integer> indeces = new TreeSet<>();
		for (int i = 0; i < domainPolynomialRing.numberOfVariables(); i++) {
			indeces.add(i + 1);
		}
		for (Polynomial<T> generator : graphIdeal.generators()) {
			Monomial monomial = generator.leadingMonomial();
			boolean indexFound = false;
			for (int i = 0; i < domainPolynomialRing.numberOfVariables(); i++) {
				if (monomial.exponents()[i] != 0) {
					indeces.remove(i + 1);
					indexFound = true;
					break;
				}
			}
			if (!indexFound) {
				continue;
			}
			for (int i = 0; i < rangePolynomialRing.numberOfVariables(); i++) {
				if (monomial.exponents()[domainPolynomialRing.numberOfVariables() + i] != 0) {
					return false;
				}
			}
		}
		return indeces.isEmpty();
	}

	public boolean isImmersion() {
		return immersionClassification().isImmersion();
	}

	public boolean isComponent() {
		return isClosedImmersion() && isOpenImmersion();
	}

	public boolean isSeparated() {
		return true;
//		AffineFiberedProduct<T> product = fiberedProduct(this, this);
//		AffineMorphism<T> diagonal = product.universalProperty(domain.identityMorphism(), domain.identityMorphism());
//		return diagonal.isClosedImmersion();
	}

	public boolean hasDenseImage() {
		AffineMorphism<T> image = closure();
		if (image.getDomain().dimension() != range.dimension()) {
			return false;
		}
		return true; // maybe (not really)
	}

	private static class AffineImmersionClassification<T extends Element<T>> {
		private boolean immersion;
		private boolean closedImmersion;
		private boolean openImmersion;
		private boolean birational;
		private CoordinateIdeal<T> closedImmersionIdeal;
		private CoordinateRingElement<T> obstruction;
		private List<TExt<T>> substitutions;
		private TranscendentalFieldExtension<T> rangeTranscendental;

		private AffineImmersionClassification() {
			this.immersion = false;
			this.closedImmersion = false;
			this.openImmersion = false;
			this.birational = false;
			this.closedImmersionIdeal = null;
			this.obstruction = null;
			this.substitutions = null;
			this.rangeTranscendental = null;
		}

		private AffineImmersionClassification(boolean closedImmersion, boolean openImmersion, boolean birational,
				CoordinateIdeal<T> closedImmersionIdeal, CoordinateRingElement<T> obstruction,
				List<TExt<T>> substitutions, TranscendentalFieldExtension<T> rangeTranscendental) {
			this.immersion = true;
			this.closedImmersion = closedImmersion;
			this.openImmersion = openImmersion;
			this.birational = birational;
			this.closedImmersionIdeal = closedImmersionIdeal;
			this.obstruction = obstruction;
			this.substitutions = substitutions;
			this.rangeTranscendental = rangeTranscendental;
		}

		public boolean isImmersion() {
			return immersion;
		}

		public boolean isClosedImmersion() {
			return closedImmersion;
		}

		public boolean isOpenImmersion() {
			return openImmersion;
		}

		public boolean isBirational() {
			return birational;
		}

		public CoordinateIdeal<T> getClosedImmersionIdeal() {
			return closedImmersionIdeal;
		}

		public CoordinateRingElement<T> getObstruction() {
			return obstruction;
		}

		public List<TExt<T>> getSubstitutions() {
			return substitutions;
		}

		public TranscendentalFieldExtension<T> getRangeTranscendental() {
			return rangeTranscendental;
		}

	}

	private AffineImmersionClassification<T> immersionClassification() {
		computeImmersionClassification();
		return immersionClassification;
	}

	private void computeImmersionClassification() {
		if (immersionClassification != null) {
			return;
		}
//		if (!range.isIntegral()) {
//			// TODO
//		}
//		AffineSimplification<T> domainSimplification = domain.simplify();
//		AffineMorphism<T> morphism = AffineMorphism.composition(domainSimplification.getInverse(), this);
		CoordinateRing<T> domainCoordinateRing = domain.getCoordinateRing();
		PolynomialRing<T> domainPolynomialRing = domainCoordinateRing.getPolynomialRing();
		CoordinateRing<T> rangeCoordinateRing = range.getCoordinateRing();
		PolynomialRing<T> rangePolynomialRing = rangeCoordinateRing.getPolynomialRing();
		PolynomialRing<T> eliminatePolynomialRing = AbstractPolynomialRing.getPolynomialRing(
				domainPolynomialRing.getRing(),
				domainPolynomialRing.numberOfVariables() + rangePolynomialRing.numberOfVariables(),
				new Monomial.EliminationOrder(Monomial.LEX, Monomial.GREVLEX,
						domainPolynomialRing.numberOfVariables()));
		AffineScheme<T> graph = new AffineScheme<>(range.getField(),
				eliminatePolynomialRing.getEmbedding(graph().getProduct().getCoordinateRing().getIdeal()).divideOut());
		CoordinateRing<T> graphCoordinateRing = graph.getCoordinateRing();
		PolynomialRing<T> graphPolynomialRing = graphCoordinateRing.getPolynomialRing();
		int[] map = new int[domainPolynomialRing.numberOfVariables() + rangePolynomialRing.numberOfVariables()];
		int[] rangeMap = new int[rangePolynomialRing.numberOfVariables()];
		TranscendentalFieldExtension<T> rangeTranscendental = new TranscendentalFieldExtension<>(range.getField(),
				rangePolynomialRing);
		TranscendentalFieldExtension<T> graphTranscendental = new TranscendentalFieldExtension<>(graph.getField(),
				graphPolynomialRing);
		List<TExt<T>> substitutions = new ArrayList<>();
		for (int i = 0; i < domainPolynomialRing.numberOfVariables(); i++) {
			map[i] = -1;
			substitutions.add(null);
		}
		for (int i = 0; i < rangePolynomialRing.numberOfVariables(); i++) {
			substitutions.add(rangeTranscendental.getVar(i + 1));
			map[domainPolynomialRing.numberOfVariables() + i] = i;
			rangeMap[i] = domainPolynomialRing.numberOfVariables() + i;
		}
		Set<Integer> variables = new TreeSet<>();
		Set<Integer> birationalVariables = new TreeSet<>();
		CoordinateIdeal<T> closedIdeal = rangeCoordinateRing.getZeroIdeal();
		graphGenerators: for (int index = graph.getCoordinateRing().getIdeal().generators().size()
				- 1; index >= 0; index--) {
			Polynomial<T> graphGenerator = graph.getCoordinateRing().getIdeal().generators().get(index);
			int variable = -1;
			for (int i = 0; i < domainPolynomialRing.numberOfVariables(); i++) {
				int degree = graphGenerator.degree(i + 1);
				if (degree > 1) {
					continue graphGenerators;
//					immersionClassification = new AffineImmersionClassification<>();
//					return;
				}
//				if (variable != -1 && i + 1 != variable && degree > 0) {
//					continue graphGenerators;
////					immersionClassification = new AffineImmersionClassification<>();
////					return;
//				}
				if (variable == -1 && degree == 1) {
					variable = i + 1;
					break;
				}
			}
			if (variable == -1) {
				closedIdeal = (CoordinateIdeal<T>) rangeCoordinateRing
						.add(rangeCoordinateRing.getIdeal(Collections.singletonList(rangeCoordinateRing
								.getEmbedding(rangePolynomialRing.getEmbedding(graphGenerator, map)))), closedIdeal);
			} else {
				UnivariatePolynomial<Polynomial<T>> collected = eliminatePolynomialRing
						.asUnivariatePolynomial(graphGenerator, variable);
				UnivariatePolynomialRing<Polynomial<T>> stackedPolynomialRing = collected.getPolynomialRing();
				CoordinateRingElement<T> linear = graphCoordinateRing
						.getEmbedding(graphPolynomialRing.fromUnivariatePolynomial(stackedPolynomialRing
								.getPolynomial(Collections.singletonList(collected.leadingCoefficient())), variable));
				CoordinateRingElement<T> constant = graphCoordinateRing
						.negative(
								graphCoordinateRing.getEmbedding(graphPolynomialRing.fromUnivariatePolynomial(
										stackedPolynomialRing.getPolynomial(
												Collections.singletonList(collected.univariateCoefficient(0))),
										variable)));
				for (int i = variable; i < domainPolynomialRing.numberOfVariables(); i++) {
					int degree = graphGenerator.degree(i + 1);
					if (degree > 0 && substitutions.get(i) == null) {
						continue graphGenerators;
					}
				}
				birationalVariables.add(variable);
				TExt<T> substitute = graphTranscendental.getElement(constant.getElement(), linear.getElement());
				substitute = rangeTranscendental.substitute(substitute, substitutions);
				if (variables.contains(variable)) {
					CoordinateRingElement<T> lhs = rangeCoordinateRing.getEmbedding(rangePolynomialRing
							.multiply(substitute.getDenominator(), substitutions.get(variable - 1).getNumerator()));
					CoordinateRingElement<T> rhs = rangeCoordinateRing.getEmbedding(rangePolynomialRing
							.multiply(substitute.getNumerator(), substitutions.get(variable - 1).getDenominator()));
					if (!lhs.equals(rhs)) {
						immersionClassification = new AffineImmersionClassification<>();
						return;
					}
				}
				CoordinateRingElement<T> inGraph = graphCoordinateRing
						.getEmbedding(graphPolynomialRing.getEmbedding(substitute.getDenominator(), rangeMap));
				boolean openVariable = false;
				if (graph.getCoordinateRing().isUnit(inGraph)) {
					variables.add(variable);
					openVariable = true;
				}
				if (openVariable || substitutions.get(variable - 1) == null) {
					substitutions.set(variable - 1, substitute);
				}
			}
		}
		boolean open = true;
		if (!closedIdeal.contains(rangeCoordinateRing.one())) {
			for (int i = 0; i < domainPolynomialRing.numberOfVariables(); i++) {
				if (!birationalVariables.contains(i + 1)) {
					immersionClassification = new AffineImmersionClassification<>();
					return;
				}
				if (!variables.contains(i + 1)) {
					open = false;
				}
			}
		}
		Polynomial<T> obstruction = rangePolynomialRing.one();
		for (int i = 0; i < domain.getCoordinateRing().getPolynomialRing().numberOfVariables(); i++) {
			obstruction = rangePolynomialRing.multiply(substitutions.get(i).getDenominator(), obstruction);
		}
		obstruction = rangePolynomialRing.radical(obstruction);
		obstruction = rangeCoordinateRing.getEmbedding(obstruction).getElement();
		obstruction = rangePolynomialRing.radical(obstruction);
		boolean closed = rangeCoordinateRing.isUnit(rangeCoordinateRing.getEmbedding(obstruction));
		open = open
				&& rangeCoordinateRing
						.intersect(closedIdeal,
								rangeCoordinateRing.getIdeal(
										Collections.singletonList(rangeCoordinateRing.getEmbedding(obstruction))))
						.equals(rangeCoordinateRing.getZeroIdeal());
		boolean birational = closedIdeal.divideOut().krullDimension() == range.dimension();
		if (birational && !closedIdeal.equals(rangeCoordinateRing.getZeroIdeal())) {
			throw new ArithmeticException("non integral range not implemented!");
		}
		immersionClassification = new AffineImmersionClassification<>(closed, open, birational, closedIdeal,
				rangeCoordinateRing.getEmbedding(obstruction),
				substitutions.subList(0, domainPolynomialRing.numberOfVariables()), rangeTranscendental);
	}

	public boolean isOpenImmersion() {
		return immersionClassification().isOpenImmersion();
	}

	public boolean isClosedImmersion() {
		return immersionClassification().isClosedImmersion();
//		if (polynomials.isEmpty()) {
//			return true;
//		}
//		CoordinateRing<T> domainCoordinateRing = domain.getCoordinateRing();
//		PolynomialRing<T> domainPolynomialRing = domainCoordinateRing.getPolynomialRing();
//		CoordinateRing<T> rangeCoordinateRing = range.getCoordinateRing();
//		PolynomialRing<T> rangePolynomialRing = rangeCoordinateRing.getPolynomialRing();
//		PolynomialRing<T> eliminatePolynomialRing = AbstractPolynomialRing.getPolynomialRing(
//				domainPolynomialRing.getRing(),
//				domainPolynomialRing.numberOfVariables() + rangePolynomialRing.numberOfVariables(),
//				new Monomial.EliminationOrder(Monomial.GREVLEX, Monomial.GREVLEX,
//						domainPolynomialRing.numberOfVariables()));
//		AffineScheme<T> graph = new AffineScheme<>(range.getField(),
//				eliminatePolynomialRing.getEmbedding(graph().getProduct().getCoordinateRing().getIdeal()).divideOut());
//		Field<T> field = domain.getField();
//		for (Polynomial<T> generator : graph.getCoordinateRing().getIdeal().generators()) {
//			boolean found = false;
//			for (int i = 0; i < domainPolynomialRing.numberOfVariables(); i++) {
//				if (generator.degree(i + 1) > 1 || (found && generator.degree(i + 1) == 1)) {
//					return false;
//				}
//				if (generator.degree(i + 1) == 1) {
//					found = true;
//				}
//			}
//		}
//		int degree = domain.degree();
//		int[] exponentLimits = new int[domainCoordinateRing.getPolynomialRing().numberOfVariables()];
//		Arrays.fill(exponentLimits, 2);
//		List<Polynomial<T>> generated = new ArrayList<>();
//		for (Polynomial<T> polynomial : polynomials) {
//			for (int i = 0; i <= degree; i++) {
//				Polynomial<T> g = domainCoordinateRing.power(domainCoordinateRing.getEmbedding(polynomial), i)
//						.getElement();
//				for (int j = 0; j < domainPolynomialRing.numberOfVariables(); j++) {
//					exponentLimits[j] = Math.max(exponentLimits[j], g.degree(j + 1) + 1);
//				}
//				generated.add(g);
//			}
//		}
//		List<Vector<T>> asVectors = new ArrayList<>();
//		for (Polynomial<T> g : generated) {
//			asVectors.add(domainPolynomialRing.asVector(g, exponentLimits));
//		}
//		Matrix<T> asMatrix = Matrix.fromColumns(asVectors);
//		for (int j = 0; j < domainPolynomialRing.numberOfVariables(); j++) {
//			if (!field.isSubModuleMember(asMatrix,
//					domainPolynomialRing.asVector(domainCoordinateRing.getVar(j + 1).getElement(), exponentLimits))) {
//				return false;
//			}
//		}
//		return true;
	}

	public CoordinateIdeal<T> closedImmersionAsIdeal() {
		if (!isClosedImmersion()) {
			throw new ArithmeticException("Not a closed immersion!");
		}
		return immersionClassification().getClosedImmersionIdeal();
	}

	public CoordinateRingElement<T> getInversionObstruction() {
		if (!isBirational()) {
			throw new ArithmeticException("Not a birational morphism!");
		}
		return immersionClassification().getObstruction();
	}

	public boolean isBirational() {
		return immersionClassification().isBirational();
	}

	public GeneralRationalFunction<T, AffinePoint<T>, AffinePoint<T>> birationalInverse() {
		if (!isBirational()) {
			throw new ArithmeticException("Not an birational function!");
		}
		CoordinateRing<T> domainCoordinateRing = domain.getCoordinateRing();
		PolynomialRing<T> domainPolynomialRing = domainCoordinateRing.getPolynomialRing();
		CoordinateRing<T> rangeCoordinateRing = range.getCoordinateRing();
		AffineMorphism<T> openImmersion = AffineMorphism
				.getOpenImmersion(range, immersionClassification().getObstruction().getElement()).getOpenImmersion();
		List<CoordinateRingElement<T>> inverse = new ArrayList<>();
		for (int i = 0; i < domainPolynomialRing.numberOfVariables(); i++) {
			TExt<T> substitute = immersionClassification().getSubstitutions().get(i);
			CoordinateRingElement<T> numerator = rangeCoordinateRing.getEmbedding(substitute.getNumerator());
			CoordinateRingElement<T> denominator = rangeCoordinateRing.getEmbedding(substitute.getDenominator());
			CoordinateRingElement<T> inducedNumerator = openImmersion.inducedMap(numerator);
			CoordinateRingElement<T> inducedDenominator = openImmersion.inducedMap(denominator);
			inverse.add(
					openImmersion.getDomain().getCoordinateRing().divideChecked(inducedNumerator, inducedDenominator));
		}
		AffineMorphism<T> inverseMorphism = new AffineMorphism<>(openImmersion.getDomain(), domain, inverse);
		return new GeneralRationalFunction<>(range, domain, inverseMorphism, 0, openImmersion, 0);
	}

	public TExt<T> birationalInverseInducedMap(CoordinateRingElement<T> t) {
		if (!isBirational()) {
			throw new ArithmeticException("Not an birational function!");
		}
		return immersionClassification().getRangeTranscendental().substitute(t.getElement(),
				immersionClassification().getSubstitutions());
	}

	public Optional<AffineMorphism<T>> inverse() {
		if (!isOpenImmersion() || !isClosedImmersion()) {
			return Optional.empty();
		}
		CoordinateRingElement<T> inverse = getInversionObstruction();
		if (!range.getCoordinateRing().isUnit(inverse)) {
			return Optional.empty();
		}
		List<Polynomial<T>> polynomials = new ArrayList<>();
		for (int i = 0; i < domain.getCoordinateRing().getPolynomialRing().numberOfVariables(); i++) {
			polynomials.add(immersionClassification().getSubstitutions().get(i).asInteger());
		}
		return Optional.of(AffineMorphism.fromPolynomials(range, domain, polynomials));
	}

	@Override
	public GeneralRationalFunction<T, AffinePoint<T>, AffinePoint<T>> restrict(AffinePoint<T> point) {
		return new GeneralRationalFunction<>(domain, range, this, 0, domain.identityMorphism(), 0);
	}

	public AffineMorphism<T> simplifyDomain() {
		AffineSimplification<T> simplification = domain.simplify();
		return composition(simplification.getInverse(), this);
	}

	public AffineMorphism<T> simplifyRange() {
		AffineSimplification<T> simplification = range.simplify();
		return composition(this, simplification.getSimplification());
	}

	public AffineMorphism<T> closure() {
		CoordinateRing<T> domainCoordinateRing = domain.getCoordinateRing();
		PolynomialRing<T> domainPolynomialRing = domainCoordinateRing.getPolynomialRing();
		CoordinateRing<T> rangeCoordinateRing = range.getCoordinateRing();
		PolynomialRing<T> rangePolynomialRing = rangeCoordinateRing.getPolynomialRing();
		// AffineFiberedProduct<T> graph = graph();
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

	public AffineFiberedProduct<T> graph() {
		return fiberedProduct(this, range.identityMorphism());
	}

	public List<AffinePoint<T>> preimageList(AffineMorphism<T> image) {
		AffineMorphism<T> preimage = preimage(image);
		if (preimage.getDomain().dimension() > 0) {
			throw new ArithmeticException("Not a finite morphism!");
		}
		List<AffinePoint<T>> result = new ArrayList<>();
		for (AffinePoint<T> preimagePoint : preimage.getDomain()) {
			result.add(preimage.evaluate(preimagePoint));
		}
		return result;
	}

	public List<AffinePoint<T>> preimageList(AffinePoint<T> image) {
		return preimageList(range.pointAsMorphism(image).restrict(SpectrumOfField.POINT).getMorphism());
	}

	/**
	 * Returns the closure of the preimage given by an injection into the range.
	 * 
	 * @return
	 */
	public AffineMorphism<T> preimage(AffineMorphism<T> image) {
		/*
		 * if (!image.isInjective()) { throw new
		 * ArithmeticException("image not an injection into range!"); }
		 */
		return fiberedProduct(this, image).projection1;
	}

	public AffineMorphism<T> preimage(AffinePoint<T> p) {
		return preimage(range.pointAsMorphism(p).restrict(SpectrumOfField.POINT).getMorphism());
	}

	public static class ExpansionResult<T extends Element<T>> {
		private AffineMorphism<T> expandedMorphism;
		private AffineMorphism<T> expandedOpenImmersion;

		private ExpansionResult(AffineMorphism<T> expandedMorphism, AffineMorphism<T> expandedOpenImmersion) {
			this.expandedMorphism = expandedMorphism;
			this.expandedOpenImmersion = expandedOpenImmersion;
		}

		public AffineMorphism<T> getExpandedMorphism() {
			return expandedMorphism;
		}

		public AffineMorphism<T> getExpandedOpenImmersion() {
			return expandedOpenImmersion;
		}

	}

	public ExpansionResult<T> expand(AffineMorphism<T> openImmersion) {
		if (!domain.equals(openImmersion.getDomain())) {
			throw new ArithmeticException("Different domain");
		}
		if (!openImmersion.isOpenImmersion()) {
			throw new ArithmeticException("Not an open immersion!");
		}
		CoordinateRing<T> domainCoordinateRing = domain.getCoordinateRing();
		CoordinateRing<T> newDomainCoordinateRing = openImmersion.getRange().getCoordinateRing();
		PolynomialRing<T> newDomainPolynomialRing = newDomainCoordinateRing.getPolynomialRing();
		List<TExt<T>> induced = new ArrayList<>();
		Polynomial<T> remainingObstruction = newDomainPolynomialRing.one();
		for (Polynomial<T> polynomial : polynomials) {
			TExt<T> inNewDomain = openImmersion
					.birationalInverseInducedMap(domainCoordinateRing.getEmbedding(polynomial));
			induced.add(inNewDomain);
			remainingObstruction = newDomainPolynomialRing.multiply(inNewDomain.getDenominator(), remainingObstruction);
		}
		AffineMorphism<T> expandedOpenImmersion = getOpenImmersion(openImmersion.getRange(), remainingObstruction)
				.getOpenImmersion();
		List<CoordinateRingElement<T>> expandedPolynomials = new ArrayList<>();
		for (TExt<T> t : induced) {
			CoordinateRingElement<T> numerator = expandedOpenImmersion
					.inducedMap(newDomainCoordinateRing.getEmbedding(t.getNumerator()));
			CoordinateRingElement<T> denominator = expandedOpenImmersion
					.inducedMap(newDomainCoordinateRing.getEmbedding(t.getDenominator()));
			expandedPolynomials
					.add(expandedOpenImmersion.getDomain().getCoordinateRing().divideChecked(numerator, denominator));
		}
		AffineMorphism<T> expandedMorphism = new AffineMorphism<>(expandedOpenImmersion.getDomain(), range,
				expandedPolynomials);
		return new ExpansionResult<>(expandedMorphism, expandedOpenImmersion);
//		AffineSimplification<T> simplifiedDomain = domain.simplify();
//		AffineMorphism<T> morphism = AffineMorphism.composition(simplifiedDomain.getInverse(), this);
//		openImmersion = AffineMorphism.composition(simplifiedDomain.getInverse(), openImmersion);
//		CoordinateRing<T> domainCoordinateRing = morphism.getDomain().getCoordinateRing();
//		PolynomialRing<T> domainPolynomialRing = domainCoordinateRing.getPolynomialRing();
//		CoordinateRing<T> rangeCoordinateRing = openImmersion.getRange().getCoordinateRing();
//		PolynomialRing<T> rangePolynomialRing = rangeCoordinateRing.getPolynomialRing();
//		PolynomialRing<T> eliminatePolynomialRing = AbstractPolynomialRing.getPolynomialRing(
//				domainPolynomialRing.getRing(),
//				domainPolynomialRing.numberOfVariables() + rangePolynomialRing.numberOfVariables(),
//				new Monomial.EliminationOrder(Monomial.GREVLEX, Monomial.GREVLEX,
//						domainPolynomialRing.numberOfVariables()));
//		Field<T> field = domain.getField();
//		AffineScheme<T> graph = new AffineScheme<>(field, eliminatePolynomialRing
//				.getEmbedding(openImmersion.graph().getProduct().getCoordinateRing().getIdeal()).divideOut());
//		List<Polynomial<T>> substitutes = new ArrayList<>();
//		for (int i = 0; i < eliminatePolynomialRing.numberOfVariables(); i++) {
//			substitutes.add(eliminatePolynomialRing.getVar(i + 1));
//		}
//		Set<Integer> variablesToRemove = new TreeSet<>();
//		for (Polynomial<T> generator : graph.getCoordinateRing().getIdeal().generators()) {
//			int variable = -1;
//			for (int i = 0; i < domainPolynomialRing.numberOfVariables(); i++) {
//				if (generator.degree(i + 1) > 1 || (variable >= 0 && generator.degree(i + 1) == 1)) {
//					throw new ArithmeticException("not an open immersion!");
//				}
//				if (generator.degree(i + 1) == 1) {
//					variable = i + 1;
//				}
//			}
//			if (variable < 0) {
//				continue;
//			}
//			boolean variableImportant = false;
//			for (Polynomial<T> polynomial : morphism.polynomials) {
//				if (polynomial.degree(variable) > 0) {
//					variableImportant = true;
//					break;
//				}
//			}
//			for (Polynomial<T> polynomial : openImmersion.polynomials) {
//				if (polynomial.degree(variable) > 0) {
//					variableImportant = true;
//					break;
//				}
//			}
//			if (!variableImportant) {
//				variablesToRemove.add(variable);
//				continue;
//			}
//			UnivariatePolynomial<Polynomial<T>> asUnivariate = eliminatePolynomialRing.asUnivariatePolynomial(generator,
//					variable);
//			int[] map = new int[eliminatePolynomialRing.numberOfVariables() - 1];
//			for (int i = 0; i < eliminatePolynomialRing.numberOfVariables(); i++) {
//				if (i + 1 < variable) {
//					map[i] = i;
//				} else if (i + 1 > variable) {
//					map[i - 1] = i;
//				}
//			}
//			Polynomial<T> leading = eliminatePolynomialRing.getEmbedding(asUnivariate.leadingCoefficient(), map);
//			Polynomial<T> constant = eliminatePolynomialRing.getEmbedding(asUnivariate.univariateCoefficient(0), map);
//			if (leading.degree() == 0) {
//				Polynomial<T> substitute = eliminatePolynomialRing
//						.divideChecked(eliminatePolynomialRing.negative(constant), leading);
//				substitutes.set(variable - 1, substitute);
//				variablesToRemove.add(variable);
//			}
//		}
//		PolynomialRing<T> expandedPolynomialRing = AbstractPolynomialRing.getPolynomialRing(field,
//				eliminatePolynomialRing.numberOfVariables() - variablesToRemove.size(), Monomial.GREVLEX);
//		int[] map = new int[eliminatePolynomialRing.numberOfVariables()];
//		int index = 0;
//		for (int i = 0; i < eliminatePolynomialRing.numberOfVariables(); i++) {
//			if (variablesToRemove.contains(i + 1)) {
//				map[i] = -1;
//				continue;
//			}
//			map[i] = index;
//			index++;
//		}
//		List<Polynomial<T>> expandedDomainGenerators = new ArrayList<>();
//		for (Polynomial<T> generator : graph.getCoordinateRing().getIdeal().generators()) {
//			Polynomial<T> substitutedGenerator = eliminatePolynomialRing.substitute(generator, substitutes);
//			boolean ignore = false;
//			for (int variable : variablesToRemove) {
//				if (substitutedGenerator.degree(variable) > 0) {
//					ignore = true;
//					break;
//				}
//			}
//			if (ignore) {
//				continue;
//			}
//			expandedDomainGenerators.add(expandedPolynomialRing.getEmbedding(substitutedGenerator, map));
//		}
//		AffineScheme<T> expandedDomain = new AffineScheme<>(domain.getField(),
//				expandedPolynomialRing.getIdeal(expandedDomainGenerators).divideOut());
//		List<Polynomial<T>> expandedPolynomials = new ArrayList<>();
//		for (Polynomial<T> polynomial : morphism.polynomials) {
//			Polynomial<T> embedded = eliminatePolynomialRing.getEmbedding(polynomial);
//			expandedPolynomials.add(expandedPolynomialRing
//					.getEmbedding(eliminatePolynomialRing.substitute(embedded, substitutes), map));
//		}
//		AffineMorphism<T> expandedMorphism = AffineMorphism.fromPolynomials(expandedDomain, morphism.getRange(),
//				expandedPolynomials);
//		List<Polynomial<T>> expandedEmbeddingPolynomials = new ArrayList<>();
//		for (Polynomial<T> polynomial : openImmersion.polynomials) {
//			Polynomial<T> embedded = eliminatePolynomialRing.getEmbedding(polynomial);
//			expandedEmbeddingPolynomials.add(expandedPolynomialRing
//					.getEmbedding(eliminatePolynomialRing.substitute(embedded, substitutes), map));
//		}
//		AffineMorphism<T> expandedImmersion = AffineMorphism.fromPolynomials(expandedDomain, openImmersion.getRange(),
//				expandedEmbeddingPolynomials);
//		AffineSimplification<T> simplification = expandedDomain.simplify();
//		expandedDomain = simplification.getSimplification().getRange();
//		expandedMorphism = AffineMorphism.composition(simplification.getInverse(), expandedMorphism);
//		expandedImmersion = AffineMorphism.composition(simplification.getInverse(), expandedImmersion);
//		return new ExpansionResult<>(expandedMorphism, expandedImmersion);
	}
}
