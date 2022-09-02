package varieties;

import fields.interfaces.Element;
import fields.interfaces.MathMap;

public interface Morphism<T extends Element<T>, S extends Element<S>, U extends Element<U>> extends MathMap<S, U> {
	Scheme<T, S> getDomain();

	Scheme<T, U> getRange();

	GeneralRationalFunction<T, S, U> restrict(S preimage);

/*
	public static class FiberedProduct<T extends Element<T>, F1 extends Element<F1>, F2 extends Element<F2>, B extends Element<B>> {
		private GluedScheme<T> product;
		private Scheme<T, F1> factor1;
		private Scheme<T, F2> factor2;
		private Scheme<T, B> base;
		private Morphism<T, GluedPoint<T>, F1> projection1;
		private Morphism<T, GluedPoint<T>, F2> projection2;
		private Morphism<T, F1, B> structural1;
		private Morphism<T, F2, B> structural2;

		public GluedScheme<T> getProduct() {
			return product;
		}

		public Scheme<T, F1> getFactor1() {
			return factor1;
		}

		public Scheme<T, F2> getFactor2() {
			return factor2;
		}

		public Scheme<T, B> getBase() {
			return base;
		}

		public Morphism<T, GluedPoint<T>, F1> getProjection1() {
			return projection1;
		}

		public Morphism<T, GluedPoint<T>, F2> getProjection2() {
			return projection2;
		}

		public Morphism<T, F1, B> getStructural1() {
			return structural1;
		}

		public Morphism<T, F2, B> getStructural2() {
			return structural2;
		}

		public<S extends Element<S>> Morphism<T, S, GluedPoint<T>> universalProperty(Morphism<T, S, F1> projection1, Morphism<T, S, F2> projection2) {
			if (!AffineMorphism.composition(projection1, structural1)
					.equals(AffineMorphism.composition(projection2, structural2))) {
				throw new ArithmeticException("Diagram does not commute!");
			}
			List<Polynomial<T>> polynomials = new ArrayList<>();
			polynomials.addAll(projection1.getPolynomials());
			polynomials.addAll(projection2.getPolynomials());
			return AffineMorphism.fromPolynomials(projection1.getDomain(), product, polynomials);
		}
	}

	public static <T extends Element<T>, F1 extends Element<F1>, F2 extends Element<F2>, B extends Element<B>> FiberedProduct<T, F1, F2, B> fiberedProduct(Morphism<T, F1, B> structural1,
			Morphism<T, F2, B> structural2) {
				FiberedProduct<T, F1, F2, B> result = new FiberedProduct<>();
		result.structural1 = structural1;
		result.structural2 = structural2;
		result.factor1 = structural1.getDomain();
		result.factor2 = structural2.getDomain();
		result.base = structural1.getRange();
		List<List<Integer>> indexTuples = new ArrayList<>();
		for (int i = 0; i < result.factor1.getAffineCover().size(); i++) {
			for (int j = 0; j < result.factor2.getAffineCover().size(); j++) {
				for (int k = 0; k < result.base.getAffineCover().size(); k++) {
					List<Integer> indeces = new ArrayList<>();
					indeces.add(i);
					indeces.add(j);
					indeces.add(k);
							indexTuples.add(indeces);
			}	}
		}
		List<AffineScheme<T>> cover = new ArrayList<>();
		List<List<AffineMorphism<T>>> atlas = new ArrayList<>();
		for (int i = 0; i < indexTuples.size(); i++) {
			List<Integer> indeces = indexTuples.get(i);
			GeneralRationalFunction<T, F1, B> restrictedStructual1 = structural1.restrict(indeces.get(0), indeces.get(2));
			GeneralRationalFunction<T, F2, B> restrictedStructual2 = structural2.restrict(indeces.get(1), indeces.get(2));
		AffineFiberedProduct<T> fiberedProduct = AffineMorphism.fiberedProduct(restrictedStructual1.getMorphism(), restrictedStructual2.getMorphism());
		}
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
*/
}
