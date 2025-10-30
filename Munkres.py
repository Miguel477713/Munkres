from typing import List, Tuple
import os

if os.name == 'nt':
    os.system('cls')  # Clear the console for Windows
elif os.name == 'posix':
    os.system('clear')  # Clear the console for Linux/Mac

Index = Tuple[int, int] #Alias stratetegy to write List[Index] insetad of Index = Tuple[int, int]

class HungarianAssignment:
    def __init__(self, costMatrix: List[List[float]], isMaximization: bool = False, showReductions: bool = False) -> None:
        """
        Initializes the solver with a square matrix.
        If isMaximization is True, converts a profit matrix to a cost matrix by inversion.
        """
        if not costMatrix or len(costMatrix) != len(costMatrix[0]):
            raise ValueError("Cost matrix must be square and non empty.")
        self.MatrixSize: int = len(costMatrix)
        self.IsMaximization: bool = isMaximization
        self.ShowReductions: bool = showReductions

        # Store original matrix for reporting
        self.OriginalMatrix: List[List[float]] = [row[:] for row in costMatrix]

        # Convert profit to cost for maximization case
        if self.IsMaximization:
            maximumValue = max(max(row) for row in costMatrix)
            self.CostMatrix: List[List[float]] = [
                [maximumValue - value for value in row] for row in costMatrix
            ]
        else:
            self.CostMatrix = [row[:] for row in costMatrix]

        if self.IsMaximization and self.ShowReductions:
            print(f"\nConverted Profit Matrix To Cost Matrix (using M = {maximumValue}):")
            self.PrintMatrix("Cost Matrix Prepared For Minimization:", self.CostMatrix)

        # Working matrix for transformations
        self.CurrentMatrix: List[List[float]] = [row[:] for row in self.CostMatrix]

        # Solution containers
        self.Assignments: List[Index] = []
        self.CoveredRows: List[bool] = [False] * self.MatrixSize
        self.CoveredColumns: List[bool] = [False] * self.MatrixSize
        self.TotalValue: float = 0.0

    # Utility printer
    def PrintMatrix(self, Title: str, Matrix: List[List[float]]) -> None:
        print(Title)
        for Row in Matrix:
            print(Row)
        print()

    def PrintCoverage(self) -> None:
        """Show which rows and columns are currently covered."""
        coveredRowList = [index + 1 for index, flag in enumerate(self.CoveredRows) if flag]
        coveredColumnList = [index + 1 for index, flag in enumerate(self.CoveredColumns) if flag]
        print(f"Covered rows: {coveredRowList}")
        print(f"Covered columns: {coveredColumnList}")

    # Step 1: Row Reduction
    def RowReduction(self) -> None:
        for rowIndex in range(self.MatrixSize):
            rowMinimum = min(self.CurrentMatrix[rowIndex])
            for columnIndex in range(self.MatrixSize):
                self.CurrentMatrix[rowIndex][columnIndex] -= rowMinimum
        if self.ShowReductions:
            self.PrintMatrix("After Row Reduction:", self.CurrentMatrix)

    # Step 2: Column Reduction
    def ColumnReduction(self) -> None:
        for columnIndex in range(self.MatrixSize):
            columnMinimum = min(self.CurrentMatrix[rowIndex][columnIndex] for rowIndex in range(self.MatrixSize))
            for rowIndex in range(self.MatrixSize):
                self.CurrentMatrix[rowIndex][columnIndex] -= columnMinimum
        if self.ShowReductions:
            self.PrintMatrix("After Column Reduction:", self.CurrentMatrix)

    def SelectIndependentZeros(self) -> List[Index]:
        """
        Greedy selection of independent zeros: pick a zero in rows with the fewest zeros first.
        This gives a valid set of independent zeros for building the minimal cover.
        """
        zeroMask: List[List[bool]] = [
            [self.CurrentMatrix[rowIndex][columnIndex] == 0 for columnIndex in range(self.MatrixSize)]
            for rowIndex in range(self.MatrixSize)
        ]
        usedRow: List[bool] = [False] * self.MatrixSize
        usedColumn: List[bool] = [False] * self.MatrixSize
        picks: List[Index] = []

        while True:
            candidateRow = -1
            candidateZeroCount = None
            for rowIndex in range(self.MatrixSize):
                if usedRow[rowIndex]:
                    continue
                count = sum(1 for columnIndex in range(self.MatrixSize)
                            if zeroMask[rowIndex][columnIndex] and not usedColumn[columnIndex])
                if count > 0 and (candidateZeroCount is None or count < candidateZeroCount):
                    candidateZeroCount = count
                    candidateRow = rowIndex
            if candidateRow == -1:
                break
            chosenColumn = None
            for columnIndex in range(self.MatrixSize):
                if zeroMask[candidateRow][columnIndex] and not usedColumn[columnIndex]:
                    chosenColumn = columnIndex
                    break
            if chosenColumn is None:
                usedRow[candidateRow] = True
                continue
            picks.append((candidateRow, chosenColumn))
            usedRow[candidateRow] = True
            usedColumn[chosenColumn] = True
        return picks

    def BuildMinimalCover(self, picks: List[Index]) -> Tuple[List[bool], List[bool]]:
        """
        Builds the minimal set of covering lines using the classic marking routine:
        1) Mark all rows that do not contain a pick
        2) While possible:
           - For each marked row, mark every column that has a zero in that row
           - For each pick located in a marked column, mark its row
        3) Covered rows are the unmarked rows, covered columns are the marked columns
        """
        zeroMask: List[List[bool]] = [
            [self.CurrentMatrix[rowIndex][columnIndex] == 0 for columnIndex in range(self.MatrixSize)]
            for rowIndex in range(self.MatrixSize)
        ]
        hasPickInRow: List[bool] = [False] * self.MatrixSize
        for rowIndex, columnIndex in picks:
            hasPickInRow[rowIndex] = True

        markedRow: List[bool] = [not hasPickInRow[rowIndex] for rowIndex in range(self.MatrixSize)]
        markedColumn: List[bool] = [False] * self.MatrixSize

        changed = True
        while changed:
            changed = False
            for rowIndex in range(self.MatrixSize):
                if not markedRow[rowIndex]:
                    continue
                for columnIndex in range(self.MatrixSize):
                    if zeroMask[rowIndex][columnIndex] and not markedColumn[columnIndex]:
                        markedColumn[columnIndex] = True
                        changed = True
            for rowIndex, columnIndex in picks:
                if markedColumn[columnIndex] and not markedRow[rowIndex]:
                    markedRow[rowIndex] = True
                    changed = True

        coveredRow: List[bool] = [not markedRow[rowIndex] for rowIndex in range(self.MatrixSize)]
        coveredColumn: List[bool] = [markedColumn[columnIndex] for columnIndex in range(self.MatrixSize)]
        return coveredRow, coveredColumn

    # Step 3: Cover Zeros With Minimal Number Of Lines
    def CoverZeros(self) -> None:
        picks = self.SelectIndependentZeros()
        coveredRow, coveredColumn = self.BuildMinimalCover(picks)
        self.Assignments = picks
        self.CoveredRows = coveredRow
        self.CoveredColumns = coveredColumn

    # Step 4: Adjust Matrix To Create More Zeros
    def AdjustMatrix(self) -> None:
        smallestUncovered = None
        uncoveredCells = []          # (rowIndex, columnIndex, value) for uncovered cells
        kPositions = []              # coordinates where value == smallestUncovered

        # Scan uncovered cells and track candidates for k
        for rowIndex in range(self.MatrixSize):
            if self.CoveredRows[rowIndex]:
                continue
            for columnIndex in range(self.MatrixSize):
                if self.CoveredColumns[columnIndex]:
                    continue
                value = self.CurrentMatrix[rowIndex][columnIndex]
                uncoveredCells.append((rowIndex, columnIndex, value))
                if smallestUncovered is None or value < smallestUncovered:
                    smallestUncovered = value

        if smallestUncovered is None:
            return

        # Identify the exact coordinates where k occurs (use 1-based in printout)
        for rowIndex, columnIndex, value in uncoveredCells:
            if value == smallestUncovered:
                kPositions.append((rowIndex + 1, columnIndex + 1))

        # Print a concise report of uncovered cells and k locations
        print("Uncovered cells (row, column, value):")
        if uncoveredCells:
            print(", ".join(f"({r+1}, {c+1}, {v})" for r, c, v in uncoveredCells))
        else:
            print("(none)")

        print(f"Adjustment applied (k = {smallestUncovered}) at positions: {kPositions}")

        # Apply the standard Hungarian adjustment
        for rowIndex in range(self.MatrixSize):
            for columnIndex in range(self.MatrixSize):
                rowIsUncovered = not self.CoveredRows[rowIndex]
                columnIsUncovered = not self.CoveredColumns[columnIndex]
                rowIsCovered = self.CoveredRows[rowIndex]
                columnIsCovered = self.CoveredColumns[columnIndex]

                if rowIsUncovered and columnIsUncovered:
                    self.CurrentMatrix[rowIndex][columnIndex] -= smallestUncovered
                elif rowIsCovered and columnIsCovered:
                    self.CurrentMatrix[rowIndex][columnIndex] += smallestUncovered

        if self.ShowReductions:
            self.PrintMatrix("After Adjustment:", self.CurrentMatrix)

    # Step 5: Repeat Coverage Until Full Cover Is Possible
    def RepeatUntilCovered(self) -> None:
        while True:
            self.CoverZeros()
            lineCount = sum(1 for flag in self.CoveredRows if flag) + sum(1 for flag in self.CoveredColumns if flag)

            #print how many lines are used each iteration
            print(f"Number of minimum lines used: {lineCount}. It shall be n={self.MatrixSize}")

            if self.ShowReductions:
                self.PrintCoverage()

            if lineCount >= self.MatrixSize:
                break

            #helpful trace before adjustment
            if self.ShowReductions:
                print("Not enough lines. Creating more zeros")

            self.AdjustMatrix()

    # Step 6: Select The Optimal Assignment
    def SelectOptimalAssignment(self) -> None:
        zeroMask: List[List[bool]] = [
            [self.CurrentMatrix[rowIndex][columnIndex] == 0 for columnIndex in range(self.MatrixSize)]
            for rowIndex in range(self.MatrixSize)
        ]
        usedRow: List[bool] = [False] * self.MatrixSize
        usedColumn: List[bool] = [False] * self.MatrixSize
        result: List[Index] = []

        # First, honor the currently selected independent zeros if they remain valid
        for rowIndex, columnIndex in self.Assignments:
            if zeroMask[rowIndex][columnIndex] and not usedRow[rowIndex] and not usedColumn[columnIndex]:
                result.append((rowIndex, columnIndex))
                usedRow[rowIndex] = True
                usedColumn[columnIndex] = True

        # Then, complete the matching by picking any available zero per free row
        for rowIndex in range(self.MatrixSize):
            if usedRow[rowIndex]:
                continue
            chosenColumn = None
            for columnIndex in range(self.MatrixSize):
                if zeroMask[rowIndex][columnIndex] and not usedColumn[columnIndex]:
                    chosenColumn = columnIndex
                    break
            if chosenColumn is not None:
                result.append((rowIndex, chosenColumn))
                usedRow[rowIndex] = True
                usedColumn[chosenColumn] = True

        # If still incomplete, fill remaining positions arbitrarily to keep one per row and column
        for rowIndex in range(self.MatrixSize):
            if usedRow[rowIndex]:
                continue
            for columnIndex in range(self.MatrixSize):
                if not usedColumn[columnIndex]:
                    result.append((rowIndex, columnIndex))
                    usedRow[rowIndex] = True
                    usedColumn[columnIndex] = True
                    break

        self.Assignments = result
        total = 0.0
        for rowIndex, columnIndex in self.Assignments:
            total += self.OriginalMatrix[rowIndex][columnIndex]
        self.TotalValue = total

    def Solve(self) -> List[List[float]]:
        """
        Runs the six steps and returns ONLY the visual assignment matrix.
        (Note: this matches current usage in RunHungarian; assignments and total are reported via Display.)
        """
        self.RowReduction()
        self.ColumnReduction()
        self.RepeatUntilCovered()
        self.SelectOptimalAssignment()

        assignmentMatrix: List[List[float]] = [[0.0 for _ in range(self.MatrixSize)] for _ in range(self.MatrixSize)]
        for rowIndex, columnIndex in self.Assignments:
            assignmentMatrix[rowIndex][columnIndex] = self.OriginalMatrix[rowIndex][columnIndex]
        return assignmentMatrix

    def Display(self) -> None:
        print("\nAssignments: ")
        for pairIndex, pairValue in enumerate(self.Assignments, start=1):
            rowIndex, columnIndex = pairValue
            print(f"  Worker {rowIndex + 1} To Task {columnIndex + 1}  Value {self.OriginalMatrix[rowIndex][columnIndex]}")
        print("\nTotal {} Value: {}".format("Maximum" if self.IsMaximization else "Minimum", self.TotalValue))


def RunHungarian() -> None:
    # demoMatrix: List[List[float]] = [
    #     [3, 5, 5, 4, 1],
    #     [2, 2, 0, 2, 2],
    #     [2, 4, 4, 1, 0],
    #     [0, 1, 1, 0, 0],
    #     [1, 2, 1, 3, 3]
    # ]
    demoMatrix: List[List[float]] = [
        [10, 9, 5],
        [9, 8, 3],
        [6, 4, 7]

    ]
    print("Task matrix:")
    for row in demoMatrix:
        print(row)
    print()

    solver = HungarianAssignment(demoMatrix, isMaximization=False, showReductions=True)
    assignMat = solver.Solve()
    solver.Display()
    print("\nAssignment Matrix With Chosen Values:")
    for row in assignMat:
        print(row)


if __name__ == "__main__":
    RunHungarian()
