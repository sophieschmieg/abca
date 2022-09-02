package varieties.curves;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.math.BigInteger;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import org.junit.jupiter.api.Test;

import fields.finitefields.FiniteField;
import fields.finitefields.FiniteField.FFE;
import fields.integers.Integers;
import fields.integers.Integers.IntE;
import fields.vectors.Matrix;
import fields.vectors.MatrixModule;
import fields.vectors.Vector;
import varieties.curves.elliptic.CompositionIsogeny;
import varieties.curves.elliptic.EllipticCurve;
import varieties.curves.elliptic.Isogeny;
import varieties.curves.elliptic.KernelIsogeny;
import varieties.curves.elliptic.KernelPointIsogeny;
import varieties.projective.ProjectivePoint;

class TorsionPointTest {

	@Test
	void testSupersingular() {
		FiniteField fq = FiniteField.getFiniteField(23 * 23);
		EllipticCurve<FFE> curve = EllipticCurve.fromJInvariant(fq, fq.getInteger(1728));
		assertEquals(BigInteger.valueOf(24 * 24), curve.getNumberOfElements());
		assertEquals(2, curve.getTorsionPointBasis(2).size());
		assertEquals(2, curve.getTorsionPointBasis(4).size());
		assertEquals(2, curve.getTorsionPointBasis(8).size());
		assertEquals(0, curve.getTorsionPointBasis(16).size());
		assertEquals(2, curve.getTorsionPointBasis(3).size());
		assertEquals(0, curve.getTorsionPointBasis(9).size());
		assertEquals(2, curve.getTorsionPointBasis(6).size());
		assertEquals(2, curve.getTorsionPointBasis(12).size());
		assertEquals(2, curve.getTorsionPointBasis(24).size());
		assertEquals(0, curve.getTorsionPointBasis(23).size());
		assertEquals(0, curve.getTorsionPointBasis(15).size());
	}

	@Test
	void testOrdinary() {
		FiniteField fq = FiniteField.getFiniteField(625);
		EllipticCurve<FFE> curve = EllipticCurve.fromJInvariant(fq, fq.getInteger(1728));
		assertEquals(BigInteger.valueOf(640), curve.getNumberOfElements());
		assertEquals(1, curve.getTorsionPointBasis(5).size());
		assertEquals(2, curve.getTorsionPointBasis(2).size());
		assertEquals(2, curve.getTorsionPointBasis(4).size());
		assertEquals(2, curve.getTorsionPointBasis(8).size());
		assertEquals(1, curve.getTorsionPointBasis(16).size());
		assertEquals(0, curve.getTorsionPointBasis(32).size());
		assertEquals(0, curve.getTorsionPointBasis(64).size());
		assertEquals(0, curve.getTorsionPointBasis(128).size());
		assertEquals(0, curve.getTorsionPointBasis(256).size());
		assertEquals(0, curve.getTorsionPointBasis(512).size());
		assertEquals(1, curve.getTorsionPointBasis(10).size());
		assertEquals(0, curve.getTorsionPointBasis(6).size());
	}

	@Test
	void testKernelIsogeny() {
		Integers z = Integers.z();
		FiniteField fq = FiniteField.getFiniteField(23 * 23);
		EllipticCurve<FFE> curve = EllipticCurve.fromJInvariant(fq, fq.getInteger(1728));
		assertEquals(BigInteger.valueOf(24 * 24), curve.getNumberOfElements());
		List<ProjectivePoint<FFE>> basis = curve.getTorsionPointBasis(24);
		Set<ProjectivePoint<FFE>> list = new TreeSet<>();
		list.add(curve.add(curve.multiply(12, basis.get(0)), curve.multiply(6, basis.get(1))));
		list.add(curve.add(curve.multiply(2, basis.get(0)), curve.multiply(6, basis.get(1))));
		list.add(curve.add(curve.multiply(3, basis.get(0)), curve.multiply(3, basis.get(1))));
		list.add(curve.add(curve.multiply(5, basis.get(0)), curve.multiply(1, basis.get(1))));
		KernelIsogeny<FFE> isogeny = new KernelIsogeny<>(curve, list);
		List<ProjectivePoint<FFE>> kernel = isogeny.kernelGenerators();
		assertEquals(2, kernel.size());
		Map<ProjectivePoint<FFE>, Vector<IntE>> expandedKernel = new TreeMap<>();
		for (int i = 0; i < 24; i++) {
			for (int j = 0; j < 24; j++) {
				ProjectivePoint<FFE> point = curve.add(curve.multiply(i, kernel.get(0)),
						curve.multiply(j, kernel.get(1)));
				expandedKernel.put(point, new Vector<>(z.getInteger(i), z.getInteger(j)));
			}
		}
		List<Vector<IntE>> asVectors = new ArrayList<>();
		for (ProjectivePoint<FFE> point : list) {
			assertTrue(expandedKernel.keySet().contains(point));
			asVectors.add(expandedKernel.get(point));
		}
		Matrix<IntE> matrix = Matrix.fromColumns(asVectors);
		MatrixModule<IntE> mm = matrix.getModule(z);
		assertTrue(mm.hasSolution(matrix, new Vector<>(z.one(), z.zero())));
		assertTrue(mm.hasSolution(matrix, new Vector<>(z.zero(), z.one())));
		Isogeny<FFE> manualConcat = curve.identity();
		EllipticCurve<FFE> domain = curve;
		for (ProjectivePoint<FFE> point : list) {
			ProjectivePoint<FFE> newKernel = manualConcat.evaluate(point);
			manualConcat = new CompositionIsogeny<>(manualConcat,
					new KernelPointIsogeny<>(domain, newKernel, domain.getOrder(newKernel).intValueExact()));
			domain = manualConcat.getRange();
		}
		assertEquals(manualConcat.getRange().jInvariant(), isogeny.getRange().jInvariant());
		for (ProjectivePoint<FFE> point : curve) {
			assertEquals(manualConcat.evaluate(point), isogeny.evaluate(point));
		}
	}
}
