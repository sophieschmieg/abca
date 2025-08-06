package fields.finitefields;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;

import fields.finitefields.Char2Field.C2E;
import fields.finitefields.FiniteField.FFE;

class Char2FieldTest {

	@Test
	void testChar2() {
		Char2Field f256 = Char2Field.fromFiniteField(FiniteField.getFiniteField(256));
		for (C2E t1 : f256) {
			for (C2E t2 : f256.getUnits()) {
				C2E div = f256.divide(t1, t2);
				assertEquals(t1, f256.multiply(div, t2));
			}
		}
	}

	@Test
	void testFiniteField() {
		FiniteField f256 = FiniteField.getFiniteField(256);
		for (FFE t1 : f256) {
			for (FFE t2 : f256.getUnits()) {
				FFE div = f256.divide(t1, t2);
				assertEquals(t1, f256.multiply(div, t2));
			}
		}
	}

}
