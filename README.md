# Final-Project-Computational-Photography
Procedure of seam carving:
- Energy Function
    - Assigning importance value to every pixel
    - Gradient magnitude
      - directional change in the intensity or color in an image
- Low-cost energy matrix
    -Cumulative minimum energy for any pixel at  (i, j) in horizontal and vertical direction
      -Considering the energies of all 3 neighboring pixels
    -Formula
      -horizontal: M(i,j) = Energy(i, j) + min(M(i-1, j-1), M(i, j-1), M(i+1, j-1))
      -vertical: M(i,j) = Energy(i, j) + min(M(i-1, j-1), M(i-1, j), M(i-1, j+1))
- Forming the low-energy seam line
    -Vertical seam line: Finding the lowest energy of each row
    -Horizontal seam line: Finding the lowest energy of each column
    -Preparing for the removing height or width based on the seam line
-Removing that seam line 
    -reducing height: Iterating each column to get rid of the horizontal seam line. 
    -reducing width: Iterating each row to get rid of the vertical seam line.

