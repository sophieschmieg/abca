package main;

public class TestField<T> {
	public TestField() {
		/*
		 * Field<Fraction<SmallIntegerRingElement>> rationalNumbers = new
		 * FieldOfFractions<SmallIntegerRingElement>(new SmallIntegerRing());
		 * PolynomialRing<Fraction<SmallIntegerRingElement>> polynomials = new
		 * PolynomialRing<Fraction<SmallIntegerRingElement>>(rationalNumbers, 3,
		 * Polynomial.GREVLEX); Map<Monomial, Fraction<SmallIntegerRingElement>> mapf =
		 * new TreeMap<Monomial, Fraction<SmallIntegerRingElement>>(Polynomial.LEX);
		 * Map<Monomial, Fraction<SmallIntegerRingElement>> mapg = new TreeMap<Monomial,
		 * Fraction<SmallIntegerRingElement>>(Polynomial.LEX); Map<Monomial,
		 * Fraction<SmallIntegerRingElement>> maph = new TreeMap<Monomial,
		 * Fraction<SmallIntegerRingElement>>(Polynomial.LEX); / *mapf.put(new
		 * Monomial(new int[] {2, 0, 0}), rationalNumbers.one()); mapf.put(new
		 * Monomial(new int[] {0, 2, 0}), rationalNumbers.one()); mapf.put(new
		 * Monomial(new int[] {0, 0, 2}), rationalNumbers.one()); mapf.put(new
		 * Monomial(new int[] {0, 0, 0}),
		 * rationalNumbers.negative(rationalNumbers.one())); mapg.put(new Monomial(new
		 * int[] {2, 0, 0}), rationalNumbers.one()); mapg.put(new Monomial(new int[] {0,
		 * 0, 2}), rationalNumbers.one()); mapg.put(new Monomial(new int[] {0, 1, 0}),
		 * rationalNumbers.negative(rationalNumbers.one())); maph.put(new Monomial(new
		 * int[] {1, 0, 0}), rationalNumbers.one()); maph.put(new Monomial(new int[] {0,
		 * 0, 1}), rationalNumbers.negative(rationalNumbers.one())); * / mapf.put(new
		 * Monomial(new int[] {2, 0, 0}), rationalNumbers.one()); mapf.put(new
		 * Monomial(new int[] {0, 2, 0}), rationalNumbers.one()); mapf.put(new
		 * Monomial(new int[] {0, 0, 2}),
		 * rationalNumbers.negative(rationalNumbers.one()));
		 * 
		 * mapg.put(new Monomial(new int[] {0, 0, 1}), rationalNumbers.one());
		 * mapg.put(new Monomial(new int[] {0, 0, 0}),
		 * rationalNumbers.negative(rationalNumbers.one()));
		 * 
		 * maph.put(new Monomial(new int[] {1, 0, 0}), rationalNumbers.one());
		 * maph.put(new Monomial(new int[] {0, 0, 1}),
		 * rationalNumbers.negative(rationalNumbers.one())); maph.put(new Monomial(new
		 * int[] {0, 0, 0}), rationalNumbers.one());
		 * Polynomial<Fraction<SmallIntegerRingElement>> f = new
		 * Polynomial<FieldOfFractions.Fraction<SmallIntegerRingElement>>(mapf,
		 * polynomials); Polynomial<Fraction<SmallIntegerRingElement>> g = new
		 * Polynomial<FieldOfFractions.Fraction<SmallIntegerRingElement>>(mapg,
		 * polynomials); Polynomial<Fraction<SmallIntegerRingElement>> h = new
		 * Polynomial<FieldOfFractions.Fraction<SmallIntegerRingElement>>(maph,
		 * polynomials); List<Polynomial<Fraction<SmallIntegerRingElement>>> basis1 =
		 * new ArrayList<Polynomial<Fraction<SmallIntegerRingElement>>>();
		 * List<Polynomial<Fraction<SmallIntegerRingElement>>> basis2 = new
		 * ArrayList<Polynomial<Fraction<SmallIntegerRingElement>>>(); basis1.add(f);
		 * basis1.add(g); basis2.add(f); basis2.add(h);
		 * PolynomialRing<Fraction<SmallIntegerRingElement>>.Ideal ideal1 =
		 * polynomials.getIdeal(basis1);
		 * PolynomialRing<Fraction<SmallIntegerRingElement>>.Ideal ideal2 =
		 * polynomials.getIdeal(basis2);
		 * PolynomialRing<Fraction<SmallIntegerRingElement>>.Ideal add =
		 * ideal1.add(ideal2); PolynomialRing<Fraction<SmallIntegerRingElement>>.Ideal
		 * mult = ideal1.multiply(ideal2);
		 * PolynomialRing<Fraction<SmallIntegerRingElement>>.Ideal inter =
		 * ideal1.intersect(ideal2);
		 * PolynomialRing<Fraction<SmallIntegerRingElement>>.Ideal sat =
		 * inter.saturate(g); System.out.println("Dimension: " + ideal1.dimension());
		 * //System.out.println("Grad: " + ideal1.degree()); System.out.println(ideal1);
		 * System.out.println("Dimension: " + ideal2.dimension());
		 * //System.out.println("Grad: " + ideal2.degree()); System.out.println(ideal2);
		 * System.out.println("Dimension: " + add.dimension());
		 * //System.out.println("Grad: " + add.degree()); System.out.println(add);
		 * System.out.println("Dimension: " + mult.dimension());
		 * //System.out.println("Grad: " + mult.degree()); System.out.println(mult);
		 * System.out.println("Dimension: " + inter.dimension());
		 * //System.out.println("Grad: " + inter.degree()); System.out.println(inter);
		 * System.out.println("Dimension: " + sat.dimension());
		 * //System.out.println("Grad: " + inter.degree()); System.out.println(sat);
		 * //Field<Fraction<SmallIntegerRingElement>> rationalNumbers = new
		 * FieldOfFractions<SmallIntegerRingElement>(new SmallIntegerRing());
		 * PolynomialRing<Fraction<SmallIntegerRingElement>> polynomials2 = new
		 * PolynomialRing<Fraction<SmallIntegerRingElement>>(rationalNumbers, 2,
		 * Polynomial.GREVLEX); Map<Monomial, Fraction<SmallIntegerRingElement>> mapell
		 * = new TreeMap<Monomial, Fraction<SmallIntegerRingElement>>(Polynomial.LEX);
		 * Map<Monomial, Fraction<SmallIntegerRingElement>> mapp = new TreeMap<Monomial,
		 * Fraction<SmallIntegerRingElement>>(Polynomial.LEX); Map<Monomial,
		 * Fraction<SmallIntegerRingElement>> mapq = new TreeMap<Monomial,
		 * Fraction<SmallIntegerRingElement>>(Polynomial.LEX);
		 * Fraction<SmallIntegerRingElement> one = rationalNumbers.one();
		 * Fraction<SmallIntegerRingElement> two = rationalNumbers.add(one, one);
		 * Fraction<SmallIntegerRingElement> three = rationalNumbers.add(two, one);
		 * 
		 * mapell.put(new Monomial(new int[] {0, 2}),
		 * rationalNumbers.negative(rationalNumbers.one())); mapell.put(new Monomial(new
		 * int[] {3, 0}), one); mapell.put(new Monomial(new int[] {2, 0}), two);
		 * mapell.put(new Monomial(new int[] {1, 0}), rationalNumbers.negative(three));
		 * 
		 * mapp.put(new Monomial(new int[] {1, 0}), two); mapp.put(new Monomial(new
		 * int[] {0, 1}), one);
		 * 
		 * mapq.put(new Monomial(new int[] {1, 0}), rationalNumbers.one()); mapq.put(new
		 * Monomial(new int[] {0, 0}), rationalNumbers.negative(three));
		 * Polynomial<Fraction<SmallIntegerRingElement>> ell = new
		 * Polynomial<FieldOfFractions.Fraction<SmallIntegerRingElement>>(mapell,
		 * polynomials2); Polynomial<Fraction<SmallIntegerRingElement>> p = new
		 * Polynomial<FieldOfFractions.Fraction<SmallIntegerRingElement>>(mapp,
		 * polynomials2); Polynomial<Fraction<SmallIntegerRingElement>> q = new
		 * Polynomial<FieldOfFractions.Fraction<SmallIntegerRingElement>>(mapq,
		 * polynomials2); List<Polynomial<Fraction<SmallIntegerRingElement>>> basis =
		 * new ArrayList<Polynomial<Fraction<SmallIntegerRingElement>>>();
		 * basis.add(ell); PolynomialRing<Fraction<SmallIntegerRingElement>>.Ideal ideal
		 * = polynomials2.getIdeal(basis);
		 * PolynomialRing<Fraction<SmallIntegerRingElement>>.Ideal satq =
		 * ideal.saturate(q); System.out.println("Dimension: " + ideal.dimension());
		 * //System.out.println("Grad: " + ideal1.degree()); System.out.println(ideal);
		 * System.out.println("Dimension: " + satq.dimension());
		 * //System.out.println("Grad: " + inter.degree()); System.out.println(satq);
		 * System.out.println(ideal.residue(p)); /
		 * *PolynomialRing<Fraction<SmallIntegerRingElement>> polynomials2 = new
		 * PolynomialRing<Fraction<SmallIntegerRingElement>>(rationalNumbers, 4,
		 * Polynomial.GREVLEX); Map<Monomial, Fraction<SmallIntegerRingElement>> mapa =
		 * new TreeMap<Monomial, Fraction<SmallIntegerRingElement>>(Polynomial.LEX);
		 * Map<Monomial, Fraction<SmallIntegerRingElement>> mapb = new TreeMap<Monomial,
		 * Fraction<SmallIntegerRingElement>>(Polynomial.LEX); Map<Monomial,
		 * Fraction<SmallIntegerRingElement>> mapc = new TreeMap<Monomial,
		 * Fraction<SmallIntegerRingElement>>(Polynomial.LEX); Map<Monomial,
		 * Fraction<SmallIntegerRingElement>> mapd = new TreeMap<Monomial,
		 * Fraction<SmallIntegerRingElement>>(Polynomial.LEX); mapa.put(new Monomial(new
		 * int[] {1, 1, 0, 0}), rationalNumbers.one()); //mapa.put(new Monomial(new
		 * int[] {0, 0, 0, 0}), rationalNumbers.negative(rationalNumbers.one()));
		 * 
		 * mapb.put(new Monomial(new int[] {0, 1, 1, 0}), rationalNumbers.one());
		 * mapb.put(new Monomial(new int[] {1, 0, 0, 1}), rationalNumbers.one());
		 * 
		 * mapc.put(new Monomial(new int[] {0, 0, 1, 0}), rationalNumbers.one());
		 * 
		 * mapd.put(new Monomial(new int[] {0, 0, 0, 1}), rationalNumbers.one());
		 * 
		 * Polynomial<Fraction<SmallIntegerRingElement>> a = new
		 * Polynomial<FieldOfFractions.Fraction<SmallIntegerRingElement>>(
		 * rationalNumbers, mapa, Polynomial.GREVLEX);
		 * Polynomial<Fraction<SmallIntegerRingElement>> b = new
		 * Polynomial<FieldOfFractions.Fraction<SmallIntegerRingElement>>(
		 * rationalNumbers, mapb, Polynomial.GREVLEX);
		 * Polynomial<Fraction<SmallIntegerRingElement>> dx = new
		 * Polynomial<FieldOfFractions.Fraction<SmallIntegerRingElement>>(
		 * rationalNumbers, mapc, Polynomial.GREVLEX);
		 * Polynomial<Fraction<SmallIntegerRingElement>> dy = new
		 * Polynomial<FieldOfFractions.Fraction<SmallIntegerRingElement>>(
		 * rationalNumbers, mapd, Polynomial.GREVLEX);
		 * List<Polynomial<Fraction<SmallIntegerRingElement>>> basis = new
		 * ArrayList<Polynomial<Fraction<SmallIntegerRingElement>>>(); basis.add(a);
		 * basis.add(b); PolynomialRing<Fraction<SmallIntegerRingElement>>.Ideal ideal =
		 * polynomials2.getIdeal(basis);
		 * PolynomialRing<Fraction<SmallIntegerRingElement>>.Ideal idealdx =
		 * polynomials2.getIdeal(Collections.singletonList(dx));
		 * PolynomialRing<Fraction<SmallIntegerRingElement>>.Ideal idealdy =
		 * polynomials2.getIdeal(Collections.singletonList(dy));
		 * System.out.println("Dimension: " + ideal.dimension());
		 * System.out.println(ideal);
		 * PolynomialRing<Fraction<SmallIntegerRingElement>>.Ideal dxnull =
		 * ideal.add(idealdx); PolynomialRing<Fraction<SmallIntegerRingElement>>.Ideal
		 * dynull = ideal.add(idealdy);
		 * PolynomialRing<Fraction<SmallIntegerRingElement>>.Ideal allnull =
		 * dxnull.add(idealdy); System.out.println(ideal.add(dxnull));
		 * System.out.println(ideal.add(dynull)); System.out.println("Dimension: " +
		 * allnull.dimension()); System.out.println(allnull);
		 */
	}
}
