import {
  ILoadReviewAcmgRules,
  ILoadReviewPublications,
} from "../../models/load-review.model";
import { IVUSSummary } from "../../models/vus-summary.model";

export interface ILoadReviewPageRequest {
  vusId: number;
}

export interface ILoadReviewPageResponse {
  variantSummary: IVUSSummary;
  publications: ILoadReviewPublications[];
  acmgRules: ILoadReviewAcmgRules[];
  classifications: string[];
}

export interface ISaveClassificationReviewRequest {
  vusId: number;
  newClassification: string;
  reason?: string;
  publicationIds?: number[];
  acmgRuleIds?: number[];
}

export interface ISaveClassificationReviewResponse {
  isSuccess: boolean;
}
