-- eq-number-format.lua
--
-- Quarto's crossref filter numbers display equations by appending
--   \qquad(N)
-- to the raw LaTeX of the equation (see eqQquad() in crossref.lua). For
-- docx output, the whole expression -- equation *and* number -- gets
-- converted to a single OMML object by Pandoc's math writer, so "(N)"
-- renders as part of the maths rather than as a normal equation-number
-- margin note.
--
-- This filter replaces each such Para with a borderless two-cell table:
-- the equation (with the \qquad(N) suffix stripped back out) centered in a
-- wide left cell, and "(N)" right-aligned in a narrow right cell --
-- mimicking the standard equation-numbering convention seen in Word
-- documents. paper-reference.docx's default "Table" style has no visible
-- borders, so this reads as plain equation + number, not as a table.
--
-- ORDERING: this only works if the filter runs *after* Quarto's own crossref
-- filter, which is what appends the \qquad(N). Quarto runs user filters first
-- by default, so paper.qmd lists the explicit `quarto` marker ahead of this
-- entry in its `filters:` block. Without it the equation arrives unnumbered
-- (with the literal "{#eq-label}" still present as a Str) and every Para
-- below falls through unchanged -- a silent no-op, not an error.
--
-- Vertical alignment of the two cells isn't controllable from here: docx cell
-- vAlign has no representation in Pandoc's AST, and the docx writer applies
-- the same "Table" style to every table, so it can't be set per-table from
-- the reference doc either. Word's default (top) applies.
--
-- Only applied to docx output; other formats use different code paths.

if FORMAT ~= "docx" then
  return {}
end

-- Fraction of the table width reserved for the number cell -- just
-- enough for "(" + two digits + ")" at the reference doc's table font size.
local NUMBER_COL_WIDTH = 0.075

local function extract_number(math_text)
  return math_text:match("\\qquad%(([^%)]+)%)%s*$")
end

local function strip_number(math_text)
  return (math_text:gsub("\\qquad%(([^%)]+)%)%s*$", ""))
end

function Para(el)
  if #el.content ~= 1 then
    return nil
  end

  -- Quarto's crossref filter wraps the labelled equation in an (otherwise
  -- empty) Span carrying the #eq-label as its identifier, so it becomes a
  -- docx bookmark that @eq-label cross-references can point at. Unwrap it
  -- to find the Math, but keep the Span around it so the bookmark --
  -- and therefore cross-references elsewhere in the document -- survive.
  local wrapper = nil
  local math_el = el.content[1]
  if math_el.t == "Span" then
    if #math_el.content ~= 1 then
      return nil
    end
    wrapper = math_el
    math_el = math_el.content[1]
  end

  if math_el.t ~= "Math" or math_el.mathtype ~= "DisplayMath" then
    return nil
  end

  local number = extract_number(math_el.text)
  if not number then
    return nil
  end
  math_el.text = strip_number(math_el.text)

  local eq_cell = pandoc.Cell({ pandoc.Plain({ wrapper or math_el }) }, pandoc.AlignCenter, 1, 1)
  local num_cell = pandoc.Cell(
    { pandoc.Plain({ pandoc.Str("(" .. number .. ")") }) },
    pandoc.AlignRight, 1, 1
  )
  local row = pandoc.Row({ eq_cell, num_cell })
  local body = {
    attr = pandoc.Attr(),
    body = { row },
    head = {},
    row_head_columns = 0,
  }

  return pandoc.Table(
    pandoc.Caption(nil),
    {
      { pandoc.AlignCenter, 1 - NUMBER_COL_WIDTH },
      { pandoc.AlignRight, NUMBER_COL_WIDTH },
    },
    pandoc.TableHead({}),
    { body },
    pandoc.TableFoot({})
  )
end
