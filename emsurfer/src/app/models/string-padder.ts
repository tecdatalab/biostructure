export class StringPadder {
    padEmdbId(strEmdbId: string): string{
        const numZeros = 4 - strEmdbId.length;
        for (let i = 0; i < numZeros; i++) { 
            strEmdbId = '0' + strEmdbId;
        }
        return strEmdbId;
    }
}
