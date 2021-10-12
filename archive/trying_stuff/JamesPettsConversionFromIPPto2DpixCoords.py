# -*- coding: utf-8 -*-
"""
James Petts' Java code to convert from the 3D coordinates in RTSTRUCTs to 2D 
in-image coordinates.

Copied from and provided by James Darcy via Slack on 04/02/2020.
"""


/**
	 * Returns the 2D coordinate within the image plane defined by the 3D patient
	 * coordinate, position, direction cosines and pixel dimensions.
	 * @param patCoord the coordinate in the patient coordinate system
	 * @param pos the image position patient
	 * @param row the direction cosines of the image row
	 * @param col the direction cosines of the image column
	 * @param pixDims the pixel dimensions
	 * @return the coordinate
	 * @throws IllegalArgumentException if any of the patient coordinate,
	 * position or direction cosines are not three element arrays
	 * @throws IllegalArgumentException if pixel dimensions is not a two element
	 * array
	 */
	public static double[] patientCoordToImageCoord(double[] patCoord, 
		double[] pos, double[] row, double[] col, double[] pixDims)
		throws IllegalArgumentException
	{
		if ((patCoord == null) || (pos == null) || (row == null) || (col == null) ||
			 (patCoord.length != 3) || (pos.length != 3) || (row.length != 3) ||
			 (col.length != 3))
		{
			throw new IllegalArgumentException(
				"Patient coordinate, image position and direction cosines must be double[3]");
		}
		if ((pixDims == null) || (pixDims.length != 2))
		{
			throw new IllegalArgumentException(
				"Pixel dimensions must be double[2]");
		}
		// Translated from James Petts' Javascript implementation
		double magRow = row[0]*row[0]+row[1]*row[1]+row[2]*row[2];
		double magCol = col[0]*col[0]+col[1]*col[1]+col[2]*col[2];
		if ((Math.abs(1.0-magRow) >= 0.01) || (Math.abs(1.0-magCol) >= 0.01))
		{
			throw new IllegalArgumentException(
				"Direction cosines must form unit vectors");
		}
		// 9 sets of simulataneous equations to choose from, choose which set to
		// solve based on the largest component of each direction cosine. This
		// avoids NaNs or floating point errors caused by dividing by very small
		// numbers and ensures a safe mapping.
		int xMaxIdx = findIndexOfMax(row);
		int yMaxIdx = findIndexOfMax(col);
		// Subtract ImagePositionPatient from coordinate
		double[] r = new double[]
		{
			patCoord[0]-pos[0], patCoord[1]-pos[1], patCoord[2]-pos[2]
		};
		// Create array to select the two simultaneous equations to solve
		double[] c = new double[]
		{
			r[xMaxIdx], col[xMaxIdx], row[xMaxIdx],
			r[yMaxIdx], row[yMaxIdx], col[yMaxIdx]
		};
		// General case: Solves the two choosen simulataneous equations to go from
		// the patient coordinate system to the image plane coordinates.
		double i = (c[0] - c[1]*c[3]/c[5]) /
			(c[2]*pixDims[0] * (1 - (c[1]*c[4])/(c[2]*c[5])));
		double j = (c[3] - c[4]*i*pixDims[0]) / (c[5]*pixDims[1]);
		return new double[] {i, j};
	}
