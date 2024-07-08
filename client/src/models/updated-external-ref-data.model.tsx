export interface IUpdatedExternalRefData {
    clinvarClassification?: string;
    clinvarClassificationLastEval?: string;
    clinvarClassificationReviewStatus?: string;
    clinvarCanonicalSpdi?: string;
    clinvarId?: number;
    clinvarVariationId?: string;
    clinvarErrorMsg?: string;
    rsid?: string;
    rsidDbsnpVerified: boolean;
    rsidDbsnpErrorMsgs: string;
    numOfPublications: number;
}