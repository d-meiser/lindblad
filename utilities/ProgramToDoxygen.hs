module Main where

import System.Environment

data LiterateSource = LiterateSource [CodeBlock]
  deriving (Show)

data CodeBlock = Comment [Line]
               | Code [Line]
  deriving (Show)

type Line = String

parseSource :: [Line] -> [CodeBlock]
parseSource [] = []
parseSource lines | isCommentStart (head lines) =
                      (parseComment lines):(parseSource $ restComment lines)
                  | otherwise =
                      (parseCode lines):(parseSource $ restCode lines)

isCommentStart :: Line -> Bool
isCommentStart ('/':'*':_) = True
isCommentStart _ = False

isCommentEnd :: Line -> Bool
isCommentEnd ('*':'/':[]) = True
isCommentEnd (_:rest) = isCommentEnd rest
isCommentEnd _ = False

dropUntil :: (a -> Bool) -> [a] -> [a]
dropUntil f = dropWhile (not . f)
takeUntil :: (a -> Bool) -> [a] -> [a]
takeUntil f = takeWhile (not . f)
parseComment :: [Line] -> CodeBlock
parseComment = Comment .
               takeUntil isCommentEnd .
               drop 1 .
               dropUntil isCommentStart
restComment :: [Line] -> [Line]
restComment = drop 1 . dropUntil isCommentEnd
parseCode :: [Line] -> CodeBlock
parseCode = Code . takeUntil isCommentStart
restCode :: [Line] -> [Line]
restCode = dropUntil isCommentStart

getCode :: String -> IO [String]
getCode s = fmap lines $ readFile s

render :: String -> String -> LiterateSource -> String
render groupName parentGroup (LiterateSource blocks) =
  "/**\n" ++
  "@defgroup " ++ groupName ++ " " ++ groupName ++ "\n" ++
  "@ingroup " ++ parentGroup ++ "\n" ++
  "@{\n" ++
  annotatedCode ++
  rawCode ++
  "@}\n" ++
  "*/\n"
  where
    annotatedCode =
      unlines $
      map renderBlock $
      drop 1 $
      filter (not . emptyBlock) blocks
    rawCode =
      unlines $
      (\ls -> "<h2>Plain source code</h2>\n\n":ls) $
      (\ls -> ["\\code\n"] ++ ls ++ ["\\endcode"]) $
      map rawContents $
      filter isCode $
      filter (not . emptyBlock) blocks
emptyBlock (Code [""]) = True
emptyBlock _ = False
renderBlock (Comment lines) = unlines lines
renderBlock (Code lines) = unlines $ ["\\code"] ++ lines ++ ["\\endcode"]
rawContents (Comment lines) = unlines lines
rawContents (Code lines) = unlines lines
isCode (Code _) = True
isCode _ = False

main :: IO ()
main =
  do
    [inputFile, groupname, parentgroup, outputFile] <- getArgs
    code <- getCode inputFile
    writeFile outputFile $
      render groupname parentgroup $ LiterateSource (parseSource code)

